
import std/[sequtils, tables, strutils, math]

type 

  Soc = float
  Soh = float
  Temperature = float
  Resistance = float
  Capacitance = float
  Voltage = float
  Current = float
  Charge = float
  Power = float
  SocTab = seq[Voltage]

  TempTab = seq[tuple[T: Temperature, f: float]]
  SohTab = seq[tuple[Q: Charge, f: float]]

  CellParam = object
    R0: Resistance
    R1: Resistance
    C1: Capacitance
    Q_bol: Charge
    socTab: SocTab
    tempTab: TempTab
    Cth: Capacitance # thermal capacitance, J/K
    Rth: Resistance # thermal resistance, K/W
    n_SOH_50: float
    charge_eff: float

  CellModel = object
    U1: Voltage

  Cell = ref object
    param: CellParam
    model: CellModel
    C: Charge # energy taken out
    Q_total: Charge
    T_ambient: Temperature
    T: Temperature
    I: Current
    U: Voltage
    U_ocv: Voltage
    U_ocv_prev: Voltage
    
  Module = object
    parallel: seq[Cell]

  Pack = object
    series: seq[Module]


proc newCell(param: CellParam): Cell =
  var cell = new Cell
  cell.param = param
  cell.model = CellModel(U1: 0.0)
  cell.C = 0.0
  cell.T = 273 + 25.0
  cell.T_ambient = 273 + 25.0
  return cell

proc Q_from_Ah(Ah: float): Charge = return Ah * 3600.0  
proc Q_to_Ah(Q: Charge): float = return Q / 3600.0


proc Q_eff(cell: Cell): Charge =
  let param = cell.param
  let T = cell.T
  let n = len(param.tempTab)
  var factor = 1.0
  if T <= param.tempTab[0].T:
    factor = param.tempTab[0].f
  elif T >= param.tempTab[n - 1].T:
    factor = param.tempTab[n - 1].f
  else:
    for i in 0 ..< n - 1:
      if param.tempTab[i].T <= T and T <= param.tempTab[i + 1].T:
        let f1 = param.tempTab[i].f
        let f2 = param.tempTab[i + 1].f
        let T1 = param.tempTab[i].T
        let T2 = param.tempTab[i + 1].T
        factor = f1 + (f2 - f1) * (T - T1) / (T2 - T1)
        break
  return param.Q_bol * factor


proc soh(cell: Cell): Soh =
  let n_SOH_50 = cell.param.n_SOH_50 # this many cycles is soh 0.5
  let cycles = cell.Q_total / cell.param.Q_bol
  return 1.0 - 0.5 * (cycles / n_SOH_50)

  
  

  


proc SOC_to_U(cp: CellParam, soc: Soc): float =
  let n = len(cp.socTab)
  if soc <= 0.0:
    return cp.socTab[0]
  if soc >= 1.0:
    return cp.socTab[n - 1]
  let idx = int(soc * (n - 1).float)
  let f1 = cp.socTab[idx]
  let f2 = cp.socTab[idx + 1]
  let soc1 = float(idx) / (n - 1).float
  let soc2 = float(idx + 1) / (n - 1).float
  return f1 + (f2 - f1) * (soc - soc1) / (soc2 - soc1)
  

proc get_soc(cell: Cell): Soc =
  if cell.C > 0.0:
    raise newException(ValueError, "Charge must be <= than zero")
  return (cell.Q_eff + cell.C) / cell.Q_eff


proc update(cell: var Cell, I: Current, dT: float): Voltage =
  let param = cell.param

  # Update the current in the cell
  cell.I = I

  # Update the equivalent circuit model and get voltage drop
  let U0 = I * param.R0
  let I_R1 = cell.model.U1 / param.R1
  let I_C1 = I - I_R1
  cell.model.U1 += dT * I_C1 / param.C1
  let dU = U0 + cell.model.U1

  # Update the cell charge
  var P_loss = 0.0
  var dC = I * dT
  if I > 0.0:
    dC *= param.charge_eff
    P_loss = cell.U * I * (1.0 - param.charge_eff)

  cell.C += dC
  cell.C = cell.C.clamp(-param.Q_bol, 0.0)
  cell.Q_total += abs(dC)

  # Update cell temperature
  let P_R0 = I * I * param.R0
  let P_R1 = I_R1 * I_R1 * param.R1
  let P_dis = P_R0 + P_R1 + P_loss
  let P_env = (cell.T - cell.T_ambient) / param.Rth
  cell.T += (P_dis - P_env) * dT / param.Cth

  # Get OCV from SOC
  let soc = cell.get_soc()
  cell.U_ocv_prev = cell.U_ocv
  cell.U_ocv = SOC_to_U(param, soc)
  cell.U = cell.U_ocv + dU
  return cell.U



let param = CellParam(
  R0: 7.0e-3,
  R1: 15.0e-3,
  C1: 1500.0,
  Cth: 25.0,
  Rth: 5.0,
  Q_bol: Q_from_Ah(2.5),
  n_SOH_50: 500,
  socTab: @[
    2.500, 2.710, 2.868, 2.972, 3.053, 3.115, 3.168, 3.212, 3.258, 3.304,
    3.347, 3.386, 3.422, 3.460, 3.484, 3.500, 3.519, 3.545, 3.571, 3.594,
    3.615, 3.638, 3.659, 3.679, 3.700, 3.722, 3.744, 3.766, 3.786, 3.805,
    3.823, 3.839, 3.855, 3.873, 3.893, 3.913, 3.932, 3.953, 3.977, 4.005,
    4.032, 4.055, 4.072, 4.081, 4.086, 4.090, 4.094, 4.100,
  ],
  tempTab: @[
    (-20.0, 0.4),
    (-10.0, 0.5),
    (  0.0, 0.7),
    ( 25.0, 1.0),
  ],
  charge_eff: 0.99,
)



proc report(cell: Cell) =
  echo $cell.I & " " & $cell.U & " " & $cell.get_soc()& " " & $(cell.T-273) & " " & $cell.soh()
  


proc discharge(cell: var Cell, I: Current) =
  if I >= 0.0:
    raise newException(ValueError, "Discharge current must be negative")
  while true:
    let U = cell.update(I, 1.0)
    cell.report()
    if U <= 2.5:
      break

proc charge(cell: var Cell, I: Current) =
  if I <= 0.0:
    raise newException(ValueError, "Charge current must be positive")
  while true:
    let U = cell.update(I, 1.0)
    cell.report()
    if U >= 4.1:
      break

var cell = newCell(param)
#cell.discharge(-5.0)

for i in 0 ..< 2:
  cell.discharge(-1.0)
  cell.charge(+1.0)



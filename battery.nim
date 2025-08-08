
import std/[sequtils, tables, strutils, math]

type 

  Soc = float
  Temperature = float
  Resistance = float
  Capacitance = float
  Voltage = float
  Current = float
  Charge = float
  Power = float
  SocTab = seq[Voltage]

  TempFactor = object
    T: Temperature
    f: float

  TempTab = seq[TempFactor]

  CellParam = object
    R0: Resistance
    R1: Resistance
    C1: Capacitance
    Q_bol: Charge
    socTab: SocTab
    tempTab: TempTab
    Cth: Capacitance # thermal capacitance, J/K
    Rth: Resistance # thermal resistance, K/W

  CellModel = object
    U1: Voltage

  Cell = ref object
    param: CellParam
    model: CellModel
    C: Charge # energy taken out
    T_ambient: Temperature
    T: Temperature
    
  Module = object
    parallel: seq[Cell]

  Pack = object
    series: seq[Module]


proc newCell(param: CellParam): Cell =
  var cell = new Cell
  cell.param = param
  cell.model = CellModel(U1: 0.0)
  cell.C = 0.0
  cell.T = 25.0
  cell.T_ambient = 25.0
  return cell

proc Q_from_Ah(Ah: float): Charge = return Ah * 3600.0  
proc Q_to_Ah(Q: Charge): float = return Q / 3600.0


proc Q_avail(cell: Cell): Charge =
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
  return param.Q_bol * factor - cell.C


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
  return (cell.param.Q_bol + cell.C) / cell.param.Q_bol


proc update(cell: var Cell, I: Current, dT: float): Voltage =
  let param = cell.param

  # Update the equivalent circuit model and get voltage drop
  let U0 = I * param.R0
  let I_R1 = cell.model.U1 / param.R1
  let I_C1 = I - I_R1
  cell.model.U1 += dT * I_C1 / param.C1
  let dU = U0 + cell.model.U1

  # Update the cell charge
  cell.C += I * dT
  cell.C = cell.C.clamp(-param.Q_bol, 0.0)
  
  # Update cell temperature
  let P_dis = I * dU
  let P_env = (cell.T - cell.T_ambient) / param.Rth
  cell.T += (P_dis - P_env) * dT / param.Cth

  # Get OCV from SOC
  let soc = cell.get_soc()
  let U_ocv = SOC_to_U(param, soc)

  return U_ocv + dU



let param = CellParam(
  R0: 3.5e-3,
  R1: 2.12e-3,
  C1: 9433,
  Cth: 25.0,
  Rth: 5.0,
  Q_bol: Q_from_Ah(2.5),
  socTab: @[
    2.710, 2.868, 2.972, 3.053, 3.115, 3.168, 3.212, 3.258, 3.304, 3.347,
    3.386, 3.422, 3.460, 3.484, 3.500, 3.519, 3.545, 3.571, 3.594, 3.615,
    3.638, 3.659, 3.679, 3.700, 3.722, 3.744, 3.766, 3.786, 3.805, 3.823,
    3.839, 3.855, 3.873, 3.893, 3.913, 3.932, 3.953, 3.977, 4.005, 4.032,
    4.055, 4.072, 4.081, 4.086, 4.090, 4.094, 4.100,
  ],
  tempTab: @[
    TempFactor(T: -20.0, f: 0.4),
    TempFactor(T: -10.0, f: 0.5),
    TempFactor(T: 0.0, f: 0.7),
    TempFactor(T: 25.0, f: 1.0),
  ]
)


var cell = newCell(param)


  
var I = 0.0

for i in 0 ..< 3600:
  I = -2.0
  if i > 60 and i < 120:
    I = -5.0
  if i > 180 and i < 240:
    I = -10.0
  let dT = 1.0
  let U = cell.update(I, dT)
  echo $I & " " & $U & " " & $cell.get_soc()& " " & $cell.T


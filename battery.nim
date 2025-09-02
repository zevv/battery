
import std/[sequtils, tables, strutils, math, times, options, os, strformat, random]

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
  Duration = float
  Interval = float
  SocTab = seq[Voltage]

  Lowpass = proc(value: float): float

  Tab[A, B] = seq[(A, B)]

  CellParam = object
    R0: Resistance
    R1: Resistance
    R2: Resistance
    C1: Capacitance
    C2: Capacitance
    Q_bol: Charge
    soc_tab: SocTab
    temp_tab: Tab[Temperature, float]
    T_R_tab: Tab[Temperature, float]
    SOH_R_tab: Tab[Soh, float]
    SOC_R_tab: Tab[Soc, float]
    Cth: Capacitance # thermal capacitance, J/K
    Rth: Resistance # thermal resistance, K/W
    n_SOH_50: float
    charge_eff: float
    peukert: float
    R_efficiency_factor: float # approximates charge efficiency drop for various chemical effects

  CellModel = object
    R0: Resistance
    R1: Resistance
    R2: Resistance
    U0: Voltage
    U1: Voltage
    U2: Voltage

  Cell = ref object
    param: CellParam
    model: CellModel
    C: Charge # energy taken out
    R: Resistance
    Q_total: Charge
    T_ambient: Temperature
    T: Temperature
    I: Current
    U: Voltage
    U_ocv: Voltage
    U_src: Voltage
    U_ocv_prev: Voltage
    soc_lowpass: Lowpass
    fd_log: File
    
  Module = object
    I_balance: Current
    parallel: seq[Cell]

  Pack = object
    series: seq[Module]

  Simulation = ref object
    pack: Pack
    time: float
    dt: float
    cycle_number: int


proc Q_from_Ah(Ah: float): Charge = return Ah * 3600.0  
proc Q_to_Ah(Q: Charge): float = return Q / 3600.0

proc interpolate[A, B](tab: Tab[A, B], x: A): B =
  let n = len(tab)
  if x <= tab[0][0]:
    return tab[0][1]
  elif x >= tab[n - 1][0]:
    return tab[n - 1][1]
  else:
    for i in 0 ..< n - 1:
      if tab[i][0] <= x and x <= tab[i + 1][0]:
        let f1 = tab[i][1]
        let f2 = tab[i + 1][1]
        let x1 = tab[i][0]
        let x2 = tab[i + 1][0]
        return f1 + (f2 - f1) * (x - x1) / (x2 - x1)
    raise newException(ValueError, "Interpolation failed")


proc mk_lowpass(alpha: float): Lowpass =
  var first = true
  var lp = 0.0
  return proc(value: float): float =
    if first:
      first = false
      lp = value
      return value
    else:
      lp = (1.0-alpha) * lp + (alpha) * value
      return lp


proc newCell(param: CellParam, idx: int): Cell =
  var cell = new Cell
  cell.param = param
  cell.model = CellModel()
  cell.C = 0.0
  cell.R = cell.param.R0 + cell.param.R1 + cell.param.R2
  cell.T = 20.0
  cell.T_ambient = 20.0
  cell.soc_lowpass = mk_lowpass(0.1)

  # Deviations

  #cell.param.Q_bol *= rand(1.00 .. 1.00)
  #cell.param.R0 *= rand(1.00 .. 1.00)
  if idx == 0:
    #cell.param.Q_bol *= 0.99
    cell.param.R0 *= 0.99
    discard

  let fname = fmt"/tmp/cell_{idx:02}.log"
  cell.fd_log = open(fname, fmWrite)
  return cell



proc get_soh(cell: Cell): Soh =
  let n_SOH_50 = cell.param.n_SOH_50 # this many cycles is soh 0.5
  let cycles = cell.Q_total / cell.param.Q_bol
  return 1.0 - 0.5 * (cycles / n_SOH_50)


proc get_Q_effective(cell: Cell): Charge =
  # Temperature factor
  let T_factor = interpolate(cell.param.temp_tab, cell.T)
  # Peukert factor
  var P_factor = 1.0
  if cell.I < 0:
    let I_ref = cell.param.Q_bol / 3600
    P_factor = pow(abs(cell.I) / I_ref, 1 - cell.param.peukert)
  return cell.param.Q_bol * T_factor * P_factor * cell.get_soh()


proc get_soc(cell: Cell): Soc =
  let Q_effective = cell.get_Q_effective()
  return (Q_effective + cell.C) / Q_effective


proc T_to_R_factor(cell: Cell): float =
  return interpolate(cell.param.T_R_tab, cell.T)


proc SOH_to_R_factor(cell: Cell): float =
  return interpolate(cell.param.SOH_R_tab, cell.get_soh())


proc SOC_to_R_factor(cell: Cell): float =
  let soc = cell.get_soc()
  return interpolate(cell.param.SOC_R_tab, soc)


proc SOC_to_U(cp: CellParam, soc: Soc): float =
  let n = len(cp.soc_tab)
  if soc <= 0.0:
    return cp.soc_tab[0]
  if soc >= 1.0:
    return cp.soc_tab[n - 1]
  let idx = int(soc * (n - 1).float)
  let f1 = cp.soc_tab[idx]
  let f2 = cp.soc_tab[idx + 1]
  let soc1 = float(idx) / (n - 1).float
  let soc2 = float(idx + 1) / (n - 1).float
  return f1 + (f2 - f1) * (soc - soc1) / (soc2 - soc1)


# Critique: This is a significant simplification. In a real cell, these factors
# affect the different sources of resistance differently:
# 
# Ohmic Resistance (R0): Primarily affected by temperature's influence on
# electrolyte and electrode conductivity.
# 
# Charge Transfer Resistance (R1): Strongly dependent on temperature
# (Arrhenius behavior) and the kinetics at the electrode-electrolyte interface,
# which varies with SOC.
# 
# Diffusion Resistance (R2): Highly dependent on SOC, as it relates to
# the concentration gradients of lithium ions within the electrodes.

proc update_R(cell: Cell) =
  let param = cell.param
  let R_factor = cell.T_to_R_factor() * cell.SOH_to_R_factor() * cell.SOC_to_R_factor()
  cell.model.R0 = param.R0 * R_factor
  cell.model.R1 = param.R1 * R_factor
  cell.model.R2 = param.R2 * R_factor
  cell.R = cell.model.R0 + cell.model.R1 + cell.model.R2


proc update(cell: Cell, I: Current, dt: Interval) =
  let param = cell.param

  # Update the current in the cell
  cell.I = I

  # Update the equivalent circuit model to calculate dU
  cell.model.U0 = I * cell.model.R0

  # Update first RC pair (charge transfer)
  let I_R1 = cell.model.U1 / cell.model.R1
  let I_C1 = I - I_R1
  cell.model.U1 += dt * I_C1 / param.C1

  # Update second RC pair (diffusion)
  let I_R2 = cell.model.U2 / cell.model.R2
  let I_C2 = I - I_R2
  cell.model.U2 += dt * I_C2 / param.C2

  let dU = cell.model.U0 + cell.model.U1 + cell.model.U2

  # Update the cell charge, taking into account the charge efficiency
  var P_loss = 0.0
  var dC = I * dt
  if I > 0.0:
    let dynamic_charge_eff = param.charge_eff - (cell.R * param.R_efficiency_factor)
    dC *= dynamic_charge_eff
    P_loss = cell.U * I * (1.0 - param.charge_eff)

  cell.C += dC
  #cell.C = cell.C.clamp(-param.Q_bol - 100, 100.0)
  cell.Q_total += abs(dC)

  # Update cell temperature
  let P_R0 = I * I * cell.model.R0
  let P_R1 = I_R1 * I_R1 * cell.model.R1
  let P_R2 = I_R2 * I_R2 * cell.model.R2
  let P_dis = P_R0 + P_R1 + P_R2 + P_loss
  let P_env = (cell.T - cell.T_ambient) / param.Rth
  cell.T += (P_dis - P_env) * dt / param.Cth

  # Get OCV from SOC
  let soc = cell.soc_lowpass(cell.get_soc())
  cell.U_ocv_prev = cell.U_ocv
  cell.U_ocv = SOC_to_U(param, soc)
  cell.U_src = cell.U_ocv + cell.model.U1 + cell.model.U2
  cell.U = cell.U_ocv + dU



let param = CellParam(
  R0:     0.03,
  R1:     0.02,
  R2:     0.01,
  C1:  3000.0,
  C2: 25000.0,
  Cth:   25.0,
  Rth:    5.0,
  Q_bol: Q_from_Ah(2.5),
  n_SOH_50: 800,
  soc_tab: @[
    2.300, 2.500, 2.710, 2.868, 2.972, 3.053, 3.115, 3.168, 3.212, 3.258,
    3.304, 3.347, 3.386, 3.422, 3.460, 3.484, 3.500, 3.519, 3.545, 3.571,
    3.594, 3.615, 3.638, 3.659, 3.679, 3.700, 3.722, 3.744, 3.766, 3.786,
    3.805, 3.823, 3.839, 3.855, 3.873, 3.893, 3.913, 3.932, 3.953, 3.977,
    4.005, 4.032, 4.055, 4.072, 4.081, 4.086, 4.090, 4.094, 4.100, 4.120,
    4.150, 4.200
  ],
  temp_tab: @[
    (-20.0, 0.4),
    (-10.0, 0.5),
    (  0.0, 0.7),
    ( 20.0, 1.0),
  ],
  T_R_tab: @[
    (-20.0, 1.5),
    (  0.0, 1.2),
    ( 20.0, 1.0)
  ],
  SOH_R_tab: @[
    (0.0, 1.5),
    (0.2, 1.3),
    (0.4, 1.2),
    (0.6, 1.1),
    (0.8, 1.05),
    (1.0, 1.0)
  ],
  SOC_R_tab: @[
    (0.0, 1.8),
    (0.1, 1.3),
    (0.2, 1.1),
    (0.4, 1.0),
    (0.6, 1.0),
    (0.8, 1.15),
    (0.9, 1.4),
    (1.0, 1.9)
  ],
  charge_eff: 0.96,
  peukert: 1.01,
  R_efficiency_factor: 5.0,
)


proc newSimulation(dt: Interval): Simulation =
  result = new Simulation
  result.dt = dt


proc report(sim: Simulation, cell: Cell, idx: int) =
  if sim.cycle_number mod 10 == 0:
    let line = fmt"{cell.I:>4.2f} {cell.U:>6.3f} {cell.get_soc():>4.3f} {cell.T:>4.2f} {cell.get_soh():>4.2f}"
    cell.fd_log.writeLine(line)


proc balance(sim: Simulation, pack: var Pack) =
  var U_min = 1e6
  for module in pack.series:
    for cell in module.parallel:
      if cell.U < U_min:
        U_min = cell.U
  for module in pack.series.mitems:
    var I_balance = 0.0
    for cell in module.parallel:
      let dU = cell.U - U_min
      if dU > 0.01:
        I_balance += -0.020
    module.I_balance = I_balance


proc run_while(sim: Simulation, I_pack: Current, condition: proc(cell: Cell): bool) =
  while true:
    var all_ok = true
    var i = 0

    sim.balance(sim.pack)

    for module in sim.pack.series:
     
      # Current in the module is pack current + balancing current
      let I_module = I_pack + module.I_balance

      var sum_U_div_R = 0.0
      var sum_1_div_R = 0.0

      # Calculate parallel resistance of all cells in the module
      var Rinv = 0.0
      for cell in module.parallel:
        cell.update_R()

        sum_U_div_R += cell.U_src / cell.model.R0
        sum_1_div_R += 1.0 / cell.model.R0

      let U_module = (I_module + sum_U_div_R) / sum_1_div_R

      for cell in module.parallel:

        let I_cell = (U_module - cell.U_src) / cell.model.R0

        cell.update(I_cell, sim.dt)
        sim.report(cell, i)
        if not condition(cell):
          all_ok = false
    if not all_ok:
          break
    sim.time += sim.dt


proc discharge(sim: Simulation, I: Current) =
  sim.run_while(I, proc(cell: Cell): bool = 
    cell.U > 2.5
  )

proc charge(sim: Simulation, I: Current) =
  sim.run_while(I, proc(cell: Cell): bool =
    cell.U < 4.2
  )

proc sleep(sim: Simulation, d: Duration) =
  let t_end = sim.time + d
  sim.run_while(0.0, proc(cell: Cell): bool =
    sim.time < t_end
  )


var sim = newSimulation(5.0)
sim.pack = Pack(
  series: @[
    Module(
      parallel: @[
        newCell(param, 0),
        newCell(param, 1),
      ]
    ),
    Module(
      parallel: @[
        newCell(param, 2),
        newCell(param, 3),
      ]
    ),
  ]
)


for i in 0 ..< 200:
  sim.cycle_number = i
  sim.discharge(-8.0)
  sim.sleep(1200)
  sim.charge(+4.0)
  sim.sleep(1200)



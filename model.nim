
import std/[sequtils, tables, strutils, times, os, strformat]
import std/[math, random, algorithm, tempfiles, complex]

type 

  Soc* = float
  Soh* = float
  Temperature* = float
  Resistance* = float
  Capacitance* = float
  Voltage* = float
  Current* = float
  Charge* = float
  Power* = float
  Duration* = float
  Interval* = float
  SocTab* = seq[Voltage]

  Tab[A, B] = seq[(A, B)]

  ArrheniusParam* = object
    A*: float # pre-exponential factor, 1/s
    Ea*: float # activation energy, J/mol

  RCParam* = object
    R*: Resistance # Ω for electrical
    C*: Capacitance # F for electrical

  RCtParam* = object
    R*: Resistance # K/W for thermal
    C*: Capacitance # J/K for thermal

  CellParam* = object
    vendor*: string
    model*: string
    RC_dc*: RCParam # DC resistance
    RC_trans*: RCParam # charge transfer
    RC_diff*: seq[RCParam] # diffusion model
    Q_bol*: Charge # nominal capacity at 1C, 3600*Ah
    I_leak_20*: Current # self-discharge current at 20°C
    soc_tab*: SocTab # OCV vs SOC
    T_cap_tab*: Tab[Temperature, float] # capacity factor vs temperature
    T_R_tab*: Tab[Temperature, float] # resistance factor vs temperature
    SOH_R_tab*: Tab[Soh, float] # resistance factor vs SOH
    SOC_R_tab*: Tab[Soc, float] # resistance factor vs SOC
    SOC_stress_Tab*: Tab[Soc, float] # 'stress' factor vs SOC
    entropy_tab*: Tab[Soc, float] # entropy coefficient vs SOC
    RCt_core*: RCtParam # RC thermal model core to case
    RCt_cell*: RCtParam # RC thermal model case to environment
    charge_eff*: float # nominal charge efficiency
    peukert*: float # Peukert exponent
    R_efficiency_factor*: float # approximates charge efficiency drop for various chemical effects
    ap_static*: ArrheniusParam # calendar aging
    ap_stress*: ArrheniusParam # cycling aging

  RCModel = object
    R: Resistance
    U: Voltage
    I_R: Current
    I_C: Current

  RCtModel = object
    T: Temperature
    P: Power

  Cell = object
    sim: Simulation
    param: CellParam
    fd_log*: File
    RC_dc: RCModel # 
    RC_trans: RCModel # charge transfer
    RC_diff: seq[RCModel] # diffusion
    Q: Charge # energy taken out
    R: Resistance
    T_env: Temperature
    RCt_core: RCtModel
    RCt_cell: RCtModel
    I_leak: Current
    I: Current
    I_lowpass: Current
    soc: Soc
    U*: Voltage # terminal voltage
    U_src: Voltage # source voltage (without R0 drop)
    soh: Soh

  Module = object
    I_balance: Current
    U: Voltage
    cells*: seq[Cell]

  Balancer = object
    I*: Current

  Pack = object
    I*: Current
    U_empty*: Voltage
    U_full*: Voltage
    modules*: seq[Module]

  Battery = object
    pack*: Pack
    balancer*: Balancer

  Simulation* = ref object
    pre_hook: proc(sim: Simulation, I_pack: Current)
    battery*: Battery
    time*: float
    time_report: int
    steps: int
    dt*: float
    cycle_number: int
    report_every_n: int


proc Q_from_Ah*(Ah: float): Charge = return Ah * 3600.0  
proc Q_to_Ah*(Q: Charge): float = return Q / 3600.0


proc interpolate[A, B](tab: Tab[A, B], x: A): B =
  let n = len(tab)
  if x <= tab[0][0]:
    return tab[0][1]
  if x >= tab[n - 1][0]:
    return tab[n - 1][1]
  for i in 0 ..< n - 1:
    if tab[i][0] <= x and x <= tab[i + 1][0]:
      let f1 = tab[i][1]
      let f2 = tab[i + 1][1]
      let x1 = tab[i][0]
      let x2 = tab[i + 1][0]
      return f1 + (f2 - f1) * (x - x1) / (x2 - x1)
  raise newException(ValueError, "Interpolation failed")


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


proc update_soc(cell: var Cell) =

  # Temperature factor
  let T_factor = interpolate(cell.param.T_cap_tab, cell.RCt_core.T)

  # Peukert factor
  var P_factor = 1.0
  let I_ref = 0.2 * cell.param.Q_bol / 3600
  if cell.I_lowpass < -I_ref:
    P_factor = pow(abs(cell.I_lowpass) / I_ref, 1 - cell.param.peukert)
  let Q_effective = cell.param.Q_bol * T_factor * P_factor * cell.soh

  cell.soc = (Q_effective + cell.Q) / Q_effective


proc update_R(cell: var Cell) =
  let param = cell.param

  let T_factor = interpolate(param.T_R_tab, cell.RCt_core.T)
  let SOH_factor = interpolate(param.SOH_R_tab, cell.soh)
  let SOC_factor = interpolate(param.SOC_R_tab, cell.soc)

  cell.RC_dc.R = param.RC_dc.R * T_factor * SOH_factor
  cell.RC_trans.R = param.RC_trans.R * T_factor * SOC_factor * SOH_factor

  for i, rc in cell.RC_diff.mpairs:
    rc.R = param.RC_diff[i].R * T_factor * SOC_factor * SOH_factor


# https://en.wikipedia.org/wiki/Arrhenius_equation
proc arrhenius(ap: ArrheniusParam, T: Temperature): float =
  return ap.A * exp(-ap.Ea / (8.314 * (T + 273.15)))


proc update_soh(cell: var Cell, dt: Interval) =
  let param = cell.param

  var soh_rate = 0.0

  # 'static' calendar aging
  soh_rate += -arrhenius(param.ap_static, cell.RCt_core.T)

  # Stress factors: power, SOC window
  let I_nominal = cell.param.Q_bol / 3600
  let power_stress = pow(abs(cell.I) / I_nominal, 1.5)
  let soc_stress = interpolate(param.SOC_stress_Tab, cell.soc)
  let degradation_rate = arrhenius(param.ap_stress, cell.RCt_core.T)
  soh_rate += -degradation_rate * power_stress * soc_stress
 
  cell.soh += soh_rate * dt

  if cell.soh <= 0.0:
    raise newException(ValueError, "Cell SOH dropped to zero")


proc update_RCt(m: var RCtModel, rp: RCtParam, P_in: Power, T_out: Temperature, dt: Interval) =
  m.P = (m.T - T_out) / rp.R
  m.T += (P_in - m.P) * dt / rp.C


proc update_temperature(cell: var Cell, I: Current, dt: Interval) =

  let param = cell.param

  # Power losses
  let P_R0 = I * I * cell.RC_dc.R
  let P_R_trans = cell.RC_trans.I_R * cell.RC_trans.I_R * cell.RC_trans.R
  let P_R_diff = cell.RC_diff.foldl(a + b.I_R * b.I_R * b.R, 0.0)
  let P_loss = cell.U * max(0, I) * (1.0 - param.charge_eff)
  let P_leak = abs(cell.I_leak) * cell.U
  let P_rev = interpolate(param.entropy_tab, cell.soc) * (cell.RCt_core.T+273.15) * I
  let P_dis = P_R0 + P_R_trans + P_R_diff + P_loss + P_leak + P_rev

  update_RCt(cell.RCt_core, param.RCt_core, P_dis, cell.RCt_cell.T, dt)
  update_RCt(cell.RCt_cell, param.RCt_cell, cell.RCt_core.P, cell.T_env, dt)

# Timestep the RC equivalent circuit model. Special case for C=0 (pure resistor).

proc update_RC(rc: var RCModel, rp: RCParam, I: Current, dt: Interval) =
  if rp.C > 0.0:
    rc.I_R = rc.U / rp.R
    rc.I_C = I - rc.I_R
    rc.U += dt * rc.I_C / rp.C
  else:
    rc.I_R = I
    rc.I_C = 0.0
    rc.U = I * rc.R


proc update_charge(cell: var Cell, I: Current, dt: Interval) =
  let param = cell.param

  # Update cell charge, taking into account the charge efficiency
  var dQ = I * dt
  if I > 0.0:
    let dynamic_charge_eff = param.charge_eff - (cell.RC_dc.R * param.R_efficiency_factor)
    dQ *= dynamic_charge_eff

  # Calculate leak current; given current is at 20°C, adjust for temperature
  cell.I_leak = param.I_leak_20 * pow(2, (cell.RCt_core.T - 20.0) / 10.0)
  dQ += cell.I_leak * dt

  cell.Q += dQ


proc update_voltage(cell: var Cell) =
  let param = cell.param
  let soc = cell.soc
  let U_ocv = SOC_to_U(param, soc)
  let U_diff = cell.RC_diff.foldl(a + b.U, 0.0)
  cell.U_src = U_ocv + cell.RC_trans.U + U_diff
  cell.U = U_ocv + cell.RC_dc.U + cell.RC_trans.U + U_diff


proc update(cell: var Cell, I: Current, dt: Interval) =
  let param = cell.param

  cell.I = I
  cell.I_lowpass = (cell.I_lowpass * 0.9) + (I * 0.1)

  update_RC(cell.RC_dc, param.RC_dc, I, dt)
  update_RC(cell.RC_trans, param.RC_trans, I, dt)
  for i, rc in cell.RC_diff.mpairs:
    update_RC(rc, param.RC_diff[i], I, dt)

  cell.update_soc()
  cell.update_voltage()
  cell.update_charge(I, dt)
  cell.update_temperature(I, dt)
  cell.update_soh(dt)


proc report(cell: var Cell) =
  let line = &"{cell.sim.time_report / 60} {cell.I:>4.2f} {cell.U:>6.3f} {cell.soc:>4.3f} {cell.RCt_core.T:>5.3f} {cell.RCt_cell.T:>5.3f} {cell.soh:>4.2f}"
  cell.fd_log.writeLine(line)

proc step*(sim: Simulation, I_pack: Current, dt: Interval): Voltage =

  if sim.pre_hook != nil:
    sim.pre_hook(sim, I_pack)

  let do_report = sim.cycle_number mod sim.report_every_n == 0

  for module in sim.battery.pack.modules.mitems:
   
    # Calculate parallel resistance of all cells in the module
    var sum_U_div_R = 0.0
    var sum_1_div_R = 0.0
    for cell in module.cells.mitems:
      cell.update_R()
      sum_U_div_R += cell.U_src / cell.RC_dc.R
      sum_1_div_R += 1.0 / cell.RC_dc.R

    # Current in the module is pack current + balancing current
    let I_module = I_pack + module.I_balance
    module.U = (I_module + sum_U_div_R) / sum_1_div_R
    result += module.U

    # Update each cell in the module with the calculated cell current
    for cell in module.cells.mitems:
      let I_cell = (module.U - cell.U_src) / cell.RC_dc.R
      cell.update(I_cell, dt)
      if do_report:
        cell.report()
        
  sim.time += dt
  inc sim.steps
  if do_report:
    sim.time_report += dt.int


proc init(cell: var Cell, sim: Simulation, param: CellParam) =
  cell.sim = sim
  cell.param = param
  cell.Q = 0.0
  cell.RCt_core.T = 20.0
  cell.RCt_cell.T = 20.0
  cell.T_env = 20.0
  cell.soh = 1.0
  cell.RC_diff = newSeq[RCModel](len(param.RC_diff))
  
  # Deviations 
  cell.param.Q_bol *= gauss(1.0, 0.01)
  cell.param.RC_dc.R *= gauss(1.0, 0.05)
  cell.param.I_leak_20 *= gauss(1.0, 0.30)
  cell.param.ap_static.A *= gauss(1.0, 0.05)
  cell.param.ap_static.Ea *= gauss(1.0, 0.05)
  cell.param.ap_stress.A *= gauss(1.0, 0.05)
  cell.param.ap_stress.Ea *= gauss(1.0, 0.05)

  # Open log file
  let fname = genTempPath("cell", "log")
  cell.fd_log = open(fname, fmReadWrite)
  removeFile(fname)
  
  # Run one step to initialize cell state
  cell.update_R()
  cell.update(0.0, 0.0)


proc init*(pack: var Pack, sim: Simulation, n_series: int, n_parallel: int, param: CellParam) =
  for i in 0 ..< n_series:
    var module = Module()
    for j in 0 ..< n_parallel:
      var cell = Cell()
      cell.init(sim, param)
      module.cells.add(cell)
    pack.modules.add(module)
  pack.U_empty = n_series.float * 2.50
  pack.U_full = n_series.float * 4.20


proc newSimulation*(dt: Interval): Simulation =
  result = new Simulation
  result.dt = dt



proc balance(sim: Simulation, I_pack: Current) =
  if I_pack > 0.0:
  
    let U_min = sim.battery.pack.modules.mapIt(it.U).min
    for module in sim.battery.pack.modules.mitems:
      let dU = module.U - U_min
      if module.U > 4.1 and dU > 0.02:
        module.I_balance = - sim.battery.balancer.I
      else:
        module.I_balance = 0.0


proc run*(sim: Simulation, fn: proc(sim: Simulation), count: int=1, n_report: int=1) =

  sim.pre_hook = balance
  sim.report_every_n = max(1, count div n_report)

  for i in 0 ..< count:
    sim.cycle_number = i
    fn(sim)
  let t = sim.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  stderr.write &"Completed {count} cycles, {sim.steps} steps, total time {days}d {hours}h {minutes}m\n"




import std/[sequtils, tables, strutils, times, options, os, strformat]
import std/[math, random, algorithm, tempfiles, complex]

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

  ArrheniusParam = object
    A: float # pre-exponential factor, 1/s
    Ea: float # activation energy, J/mol

  RCParam = object
    R: Resistance # Ω for electrical
    C: Capacitance # F for electrical

  RCtParam = object
    R: Resistance # K/W for thermal
    C: Capacitance # J/K for thermal

  CellParam = object
    vendor: string
    model: string
    RC_dc: RCParam # DC resistance
    RC_trans: RCParam # charge transfer
    RC_diff: seq[RCParam] # diffusion model
    Q_bol: Charge # nominal capacity at 1C, 3600*Ah
    I_leak_20: Current # self-discharge current at 20°C
    soc_tab: SocTab # OCV vs SOC
    temp_tab: Tab[Temperature, float] # capacity factor vs temperature
    T_R_tab: Tab[Temperature, float] # resistance factor vs temperature
    SOH_R_tab: Tab[Soh, float] # resistance factor vs SOH
    SOC_R_tab: Tab[Soc, float] # resistance factor vs SOC
    SOC_stress_Tab: Tab[Soc, float] # 'stress' factor vs SOC
    entropy_tab: Tab[Soc, float] # entropy coefficient vs SOC
    RCcore: RCtParam # RC thermal model core to case
    RCcase: RCtParam # RC thermal model case to environment
    charge_eff: float # nominal charge efficiency
    peukert: float # Peukert exponent
    R_efficiency_factor: float # approximates charge efficiency drop for various chemical effects
    ap_static: ArrheniusParam # calendar aging
    ap_stress: ArrheniusParam # cycling aging

  RCModel = object
    R: Resistance
    U: Voltage
    I_R: Current
    I_C: Current

  Cell = object
    id: int
    sim: Simulation
    param: CellParam
    RC_dc: RCModel # 
    RC_trans: RCModel # charge transfer
    RC_diff: seq[RCModel] # diffusion
    Q: Charge # energy taken out
    Q_total: Charge
    R: Resistance
    T_env: Temperature
    T_core: Temperature
    T_case: Temperature
    I: Current
    I_lowpass: Current
    U: Voltage # terminal voltage
    I_leak: Current
    U_src: Voltage # source voltage (without R0 drop)
    fd_log: File
    soh: Soh
    dsoh: float

  Module = object
    I_balance: Current
    U: Voltage
    cells: seq[Cell]

  Balancer = object
    I: Current

  Pack = object
    U_empty: Voltage
    U_full: Voltage
    modules: seq[Module]
    balancer: Balancer

  Simulation = ref object
    pack: Pack
    time: float
    time_report: int
    dt: float
    cycle_number: int
    report_every_n: int


proc Q_from_Ah(Ah: float): Charge = return Ah * 3600.0  
proc Q_to_Ah(Q: Charge): float = return Q / 3600.0


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


proc get_Q_effective(cell: Cell): Charge =
  # Temperature factor
  let T_factor = interpolate(cell.param.temp_tab, cell.T_core)
  # Peukert factor
  var P_factor = 1.0
  let I_ref = 0.2 * cell.param.Q_bol / 3600
  if cell.I_lowpass < -I_ref:
    P_factor = pow(abs(cell.I_lowpass) / I_ref, 1 - cell.param.peukert)
  return cell.param.Q_bol * T_factor * P_factor * cell.soh


proc get_soc(cell: Cell): Soc =
  let Q_effective = cell.get_Q_effective()
  return (Q_effective + cell.Q) / Q_effective


proc T_to_R_factor(cell: Cell): float =
  return interpolate(cell.param.T_R_tab, cell.T_core)


proc SOH_to_R_factor(cell: Cell): float =
  return interpolate(cell.param.SOH_R_tab, cell.soh)


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


proc update_R(cell: var Cell) =
  let param = cell.param

  let T_factor = cell.T_to_R_factor()
  let SOH_factor = cell.SOH_to_R_factor()
  let SOC_factor = cell.SOC_to_R_factor()

  cell.RC_dc.R = param.RC_dc.R * T_factor * SOH_factor
  cell.RC_trans.R = param.RC_trans.R * T_factor * SOC_factor * SOH_factor

  for i, rc in cell.RC_diff.mpairs:
    rc.R = param.RC_diff[i].R * T_factor * SOC_factor * SOH_factor


proc update_soh(cell: var Cell, dt: Interval) =
  let param = cell.param

  # Arrhenius equation to calculate a degradation factor based on the cell's
  # current temperature
  proc arrhenius(ap: ArrheniusParam, T: Temperature): float =
    return ap.A * exp(-ap.Ea / (8.314 * (T + 273.15)))

  cell.dsoh = 0.0

  # 'static' calendar aging
  cell.dsoh += -arrhenius(param.ap_static, cell.T_core)

  # 'stress' aging due to cycling
  let dC = abs(cell.I) * dt
  if dC > 0:
    # Stress factors: power, SOC window
    let I_nominal = cell.param.Q_bol / 3600
    let power_stress = pow(abs(cell.I) / I_nominal, 1.5)
    let soc_stress = interpolate(param.SOC_stress_Tab, cell.get_soc())
    let degradation_rate = arrhenius(param.ap_stress, cell.T_core)
    cell.dsoh += -degradation_rate * power_stress * soc_stress
  
  cell.soh += cell.dsoh * dt

  if cell.soh <= 0.0:
    raise newException(ValueError, "Cell SOH dropped to zero")


proc update_temperature(cell: var Cell, I: Current, dt: Interval) =

  let param = cell.param
  let P_R0 = I * I * cell.RC_dc.R
  let P_R_trans = cell.RC_trans.I_R * cell.RC_trans.I_R * cell.RC_trans.R
  let P_R_diff = cell.RC_diff.foldl(a + b.I_R * b.I_R * b.R, 0.0)
  let P_loss = cell.U * max(0, I) * (1.0 - param.charge_eff)
  let P_leak = abs(cell.I_leak) * cell.U
  let P_rev = interpolate(param.entropy_tab, cell.get_soc()) * (cell.T_core+273.15) * I
  let P_dis = P_R0 + P_R_trans + P_R_diff + P_loss + P_leak + P_rev

  let P_core_to_case = (cell.T_core - cell.T_case) / param.RCcore.R
  cell.T_core += (P_dis - P_core_to_case) * dt / param.RCcore.C
  let P_case_to_env = (cell.T_case - cell.T_env) / param.RCcase.R
  cell.T_case += (P_core_to_case - P_case_to_env) * dt / param.RCcase.C


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
  # Update the cell charge, taking into account the charge efficiency
  var dC = I * dt
  if I > 0.0:
    let dynamic_charge_eff = param.charge_eff - (cell.RC_dc.R * param.R_efficiency_factor)
    dC *= dynamic_charge_eff

  # Calculate leak current; given current is at 20°C, adjust for temperature
  cell.I_leak = param.I_leak_20 * pow(2, (cell.T_core - 20.0) / 10.0)
  dC += cell.I_leak * dt

  cell.Q += dC
  cell.Q_total += abs(dC)


proc update_voltage(cell: var Cell) =
  let param = cell.param
  let soc = cell.get_soc()
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

  cell.update_voltage()
  cell.update_charge(I, dt)
  cell.update_temperature(I, dt)
  cell.update_soh(dt)


var idx = 0

proc newCell(sim: Simulation, param: CellParam): Cell =
  var cell = Cell()
  cell.sim = sim
  cell.param = param
  cell.Q = 0.0
  cell.T_core = 20.0
  cell.T_case = 20.0
  cell.T_env = 20.0
  cell.soh = 1.0
  cell.RC_diff = newSeq[RCModel](len(param.RC_diff))

  let fname = genTempPath("cell", "log")
  cell.fd_log = open(fname, fmReadWrite)
  removeFile(fname)
  inc idx

  # Run one step to initialize cell state
  cell.update_R()
  cell.update(0.0, 0.0)

  # Deviations
  cell.param.Q_bol *= rand(0.9 .. 1.00)
  cell.param.RC_dc.R *= rand(0.95 .. 1.00)

  return cell


proc newSimulation(dt: Interval): Simulation =
  result = new Simulation
  result.dt = dt


proc gen_gnuplot(sim: Simulation, fname: string) =
  let fd = open(fname, fmWrite)
  proc l(s: string) =
    fd.write(s & "\n")

  l("#!/usr/bin/gnuplot -p")
  l("")
  l("reset")
  l("set grid")
  l("set key off")
  l("set multiplot layout 5, 1")
  l("set lmargin at screen 0.08")
  #l("set noxtics")
  l("# 0 margin between multiplots")
  l("set tmargin 1")
  l("set bmargin 1")
  l("set offsets graph 0, 0, 0.05, 0.05")
  l("")

  # Emit inline data blocks for all cells

  for mid, module in sim.pack.modules:
    for cid, cell in module.cells:
      l(&"$cell_{mid:02}_{cid:02} << EOD")
      cell.fd_log.setFilePos(0)
      for line in cell.fd_log.lines:
        l(line)
      l("EOD")
      l("")

  # Emit plotting commands

  proc gen_graph(gs: openArray[tuple[col: int, ylabel: string]], pres: openArray[string]) =
    var ts: seq[string]
    for i, g in gs:
      let lt = if i == 0: "1" else: "2"
      for mid, module in sim.pack.modules:
        for cid, cell in module.cells:
          ts.add(&""" $cell_{mid:02}_{cid:02} u 1:{g.col} w l dt {lt} """ )
    for pre in pres:
      l(pre)
    l(&"""set ylabel "{gs[0].ylabel}"""")
    l(&"""plot {ts.join(", ")}""")
    l("")

  gen_graph([ (2, "I (A)",      )], [ "unset yrange" ])
  gen_graph([ (3, "U (V)",      )], [ "set yrange [2.3:4.4]" ])
  gen_graph([ (4, "SOC (%)",    )], [ "set yrange [-0.1:1.1]" ])
  gen_graph([ (5, "T_core (°C)",),
              (6, "T_case (°C)",)], [ "set yrange [19:30] " ])
  gen_graph([ (7, "SOH (%)",    )], [ "set yrange [-0.1:1.1]" ])

  fd.write("unset multiplot\n")
  fd.close()


proc report(cell: var Cell) =
  let line = &"{cell.sim.time_report} {cell.I:>4.2f} {cell.U:>6.3f} {cell.get_soc():>4.3f} {cell.T_core:>5.3f} {cell.T_case:>5.3f} {cell.soh:>4.2f} {log10(-cell.dsoh):>g}"
  cell.fd_log.writeLine(line)
  
    # bitline
    #cell.fd_log.writeLine(&"{cell.sim.time:f} g cell{cell.id}.I {cell.I:f}")
    #cell.fd_log.writeLine(&"{cell.sim.time:f} g cell{cell.id}.U {cell.U:f}")
    #cell.fd_log.writeLine(&"{cell.sim.time:f} g cell{cell.id}.SOC {cell.get_soc():f}")
    #cell.fd_log.writeLine(&"{cell.sim.time:f} g cell{cell.id}.T {cell.T_core:f}")
    #cell.fd_log.writeLine(&"{cell.sim.time:f} g cell{cell.id}.SOH {cell.soh:f}")



proc balance(sim: Simulation, pack: var Pack) =
  let U_min = pack.modules.mapIt(it.U).min
  for module in pack.modules.mitems:
    let dU = module.U - U_min
    if module.U > 3.6 and dU > 0.01:
      module.I_balance = - pack.balancer.I
    else:
      module.I_balance = 0.0


proc step(sim: Simulation, I_pack: Current, dt: Interval): Voltage =

  sim.balance(sim.pack)
  
  let do_report = sim.cycle_number mod sim.report_every_n == 0

  for module in sim.pack.modules.mitems:
   
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
  if do_report:
    sim.time_report += dt.int


proc newPack(sim: Simulation, n_series: int, n_parallel: int, param: CellParam): Pack =
  result = Pack()
  for i in 0 ..< n_series:
    var module = Module()
    for j in 0 ..< n_parallel:
      module.cells.add(sim.newCell(param))
    result.modules.add(module)
  result.U_empty = n_series.float * 2.50
  result.U_full = n_series.float * 4.20


proc discharge(sim: Simulation, I: Current, U_min: Voltage) =
  while true:
    var U = sim.step(I, sim.dt)
    if U < U_min:
      return
    for module in sim.pack.modules:
      for cell in module.cells:
        if cell.U < 2.5:
          return

proc charge(sim: Simulation, I: Current, U_max: Voltage) =
  while true:
    var U_pack = sim.step(I, sim.dt)
    if U_pack > U_max:
      break


proc charge_CC_CV(sim: Simulation, I_set: Current, U_set: Voltage) =
  
  # PID controller constants for voltage regulation
  let kP = 0.3
  let kI = 2.5

  let t_max = sim.time + 6 * 3600
  var I_pack = I_set
  var err_int = 0.0

  while I_pack > I_set * 0.05:
    let I = clamp(I_pack, 0.0, I_set)
    var U_pack = sim.step(I, sim.dt)

    let err = U_set - U_pack
    err_int += err * sim.dt
    err_int = clamp(err_int, -2.0, 2.0)
  
    I_pack = (kP * err) + (kI * err_int)

    if sim.time > t_max:
      raise newException(ValueError, "CC/CV charge timeout")
    

proc sleep(sim: Simulation, d: Duration) =
  let t_end = sim.time + d
  while sim.time < t_end:
    discard sim.step(0.0, sim.dt)



let param = CellParam(
  vendor: "Samsung",
  model:  "INR18650-32E",
  RC_dc:      RCParam(R: 0.025),
  RC_trans: RCParam(R: 0.015, C:  4000.0),
  RC_diff: @[ 
    RCParam(R: 0.030, C:  30_000),
    RCParam(R: 0.015, C:  90_000),
    RCParam(R: 0.008, C: 270_000),
    RCParam(R: 0.007, C: 810_000),
    RCParam(R: 0.006, C:1600_000), # τ ≈ 3240s (~54 min)
  ],
  RCcore: RCtParam(R: 2.500, C:   150.0),
  RCcase: RCtParam(R: 5.000, C:    30.0),
  Q_bol: Q_from_Ah(3.2),
  I_leak_20: -0.14e-3,
  soc_tab: @[
    2.300, 2.500, 2.710, 2.868, 2.972, 3.053, 3.115, 3.168, 3.212, 3.258,
    3.304, 3.347, 3.386, 3.422, 3.460, 3.484, 3.500, 3.519, 3.545, 3.571,
    3.594, 3.615, 3.638, 3.659, 3.679, 3.700, 3.722, 3.744, 3.766, 3.786,
    3.805, 3.823, 3.839, 3.855, 3.873, 3.893, 3.913, 3.932, 3.953, 3.977,
    4.005, 4.032, 4.055, 4.072, 4.081, 4.086, 4.090, 4.094, 4.100, 4.120,
    4.150, 4.200, 4.250
  ],
  temp_tab: @[
    (-20.0, 0.60),
    (-10.0, 0.75),
    (  0.0, 0.88),
    ( 25.0, 1.00),
    ( 40.0, 1.02)
  ],
  T_R_tab: @[
    (-20.0, 3.0),
    (  0.0, 1.8),
    ( 25.0, 1.0)
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
  SOC_stress_Tab: @[
    (0.0, 3.0),
    (0.1, 1.5),
    (0.2, 1.0),
    (0.4, 1.0),
    (0.6, 1.0),
    (0.8, 1.0),
    (0.9, 1.5),
    (1.0, 3.0)
  ],
  # file:///home/ico/Downloads/energies-12-02685.pdf, figure 5
  entropy_tab: @[
    (0.0,  0.0002),
    (0.2,  0.0001),
    (0.5, -0.0001),
    (0.8, -0.0003),
    (1.0,  0.0004)
  ],
  charge_eff: 0.96,
  peukert: 1.03,
  R_efficiency_factor: 5.0,
  ap_static: ArrheniusParam(
    A: 5.0,
    Ea: 55.0e3,
  ),
  ap_stress: ArrheniusParam(
    A: 300,
    Ea: 54.0e3,
  )
)



proc test_cycle(sim: Simulation) =
  sim.sleep(600)
  sim.discharge(-8.0, sim.pack.U_empty)
  sim.sleep(3600)
  #sim.charge(+4.0, sim.pack.U_full)
  sim.charge_CC_CV(+4.0, sim.pack.U_full)
  sim.sleep(600)


proc test_sleep(sim: Simulation) =
  sim.sleep(24 * 3600)


proc test_commute(sim: Simulation) =
  proc drive(sim: Simulation, I: Current, d: Duration) =
    let t_end = sim.time + d
    while sim.time < t_end:
      discard sim.step(I, sim.dt)
  proc commute() =
    var commute_time_remaining = 40 * 60.0
    for i in 1..5:
      let drive_t = rand(45.0 .. 90.0)
      let drive_I = rand(-5.0 .. -3.0)
      sim.drive(drive_I, drive_t)
      commute_time_remaining -= drive_t
      let stop_t = rand(15.0 .. 45.0)
      sim.drive(rand(0.2 .. 0.5), stop_t)
      commute_time_remaining -= stop_t
    let highway_t = 15.0 * 60.0
    sim.drive(-8.0, highway_t)
    commute_time_remaining -= highway_t
    sim.drive(-4.0, 2 * 60)
    commute_time_remaining -= 2 * 60
    sim.drive(2.0, 30)
    commute_time_remaining -= 30
    sim.drive(0.0, commute_time_remaining)
  sim.sleep(6 * 3600)
  echo "--- Morning Commute ---"
  commute()
  echo "--- Parking at Work (8h 20m) ---"
  sim.sleep(8 * 3600 + 20 * 60)
  echo "--- Afternoon Commute ---"
  commute()
  echo "--- Parking at Home (1h 20m) ---"
  sim.sleep(1 * 3600 + 20 * 60)
  echo "--- Evening Charge ---"
  sim.charge_CC_CV(+4.0, sim.pack.U_full)
  echo "--- Overnight Rest ---"
  let total_cycle_duration = 24.0 * 3600.0
  let time_into_this_cycle = sim.time mod total_cycle_duration
  let remaining_time = total_cycle_duration - time_into_this_cycle
  if remaining_time > 0:
    sim.sleep(remaining_time)

proc test_EIS_f(sim: Simulation, freq: float) =
  let cycles = 3
  let steps_per_cycle = 20
  let steps = cycles * steps_per_cycle
  let dt = 1.0 / (freq * steps_per_cycle.float)
  stderr.write(&"EIS {freq:>8.3f} Hz, dt={dt:>6.4g} s, steps={steps}\n")
  let I_amp = 0.010
  var Zr = 0.0
  var Zi = 0.0
  var n = 0
  let t_start = sim.time
  while n < steps:
    let t = sim.time - t_start
    let ref_sin = sin(TAU * freq * t)
    let ref_cos = cos(TAU * freq * t)
    var U = sim.step(ref_sin * I_amp, dt)
    Zr  += U * ref_sin / steps.float
    Zi += U * ref_cos / steps.float
    inc n
  let pha = arctan2(Zi, Zr) * 180.0 / PI
  let mag = sqrt(Zr*Zr + Zi*Zi)
  echo freq, " ", Zr, " ", Zi, " ", mag, " ", pha


proc test_EIS_f2(sim: Simulation, freq: float) =
  let w = TAU * freq
  var Z = complex(param.RC_dc.R)
  Z += param.RC_trans.R / complex(1.0, w * param.RC_trans.R * param.RC_trans.C)
  for rc in param.RC_diff:
    Z += rc.R / complex(1.0, w * rc.R * rc.C)
  let mag = abs(Z)
  let pha = arctan2(Z.im, Z.re) * 180.0 / PI
  echo freq, " ", Z.re, " ", Z.im, " ", mag, " ", pha


proc test_EIS(sim: Simulation) =
  var f = 0.0006
  while f < 100:
    sim.charge_CC_CV(+4.0, sim.pack.U_full)
    sim.sleep(3600)
    sim.test_EIS_f(f)
    f *= 1.5

proc run(sim: Simulation, fn: proc(sim: Simulation), count: int=1, n_report: int=1) =
  sim.report_every_n = max(1, count div n_report)
  for i in 0 ..< count:
    sim.cycle_number = i
    fn(sim)
  let t = sim.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  stderr.write &"Completed {count} cycles, total time {days}d {hours}h {minutes}m\n"


var sim = newSimulation(10.0)
sim.pack = sim.newPack(n_series=2, n_parallel=2, param)
sim.pack.balancer.I = 0.200

sim.run(test_cycle, count=300, n_report=3)
#sim.run(test_EIS)
sim.gen_gnuplot("battery.gp")


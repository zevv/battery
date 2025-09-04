
import std/[sequtils, tables, strutils, math, times, options, os, strformat, random, algorithm]

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

  CellParam = object
    vendor: string
    model: string
    R0: Resistance
    R1: Resistance
    R2: Resistance
    C1: Capacitance
    C2: Capacitance
    Q_bol: Charge
    I_leak_20: Current
    soc_tab: SocTab
    temp_tab: Tab[Temperature, float]
    T_R_tab: Tab[Temperature, float]
    SOH_R_tab: Tab[Soh, float]
    SOC_R_tab: Tab[Soc, float]
    SOC_stress_Tab: Tab[Soc, float]
    Cth: Capacitance # thermal capacitance, J/K
    Rth: Resistance # thermal resistance, K/W
    charge_eff: float
    peukert: float
    R_efficiency_factor: float # approximates charge efficiency drop for various chemical effects
    ap_static: ArrheniusParam # calendar aging
    ap_stress: ArrheniusParam # cycling aging


  CellModel = object
    R0: Resistance
    R1: Resistance
    R2: Resistance
    U0: Voltage
    U1: Voltage
    U2: Voltage

  Cell = ref object
    id: int
    sim: Simulation
    param: CellParam
    model: CellModel
    C: Charge # energy taken out
    R: Resistance
    Q_total: Charge
    T_ambient: Temperature
    T: Temperature
    I: Current
    I_lowpass: Current
    U: Voltage # terminal voltage
    U_src: Voltage # source voltage (without R0 drop)
    fd_log: File
    soh: Soh
    dsoh: float

  Module = object
    I_balance: Current
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
    cells: seq[Cell]
    time: float
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
  let T_factor = interpolate(cell.param.temp_tab, cell.T)
  # Peukert factor
  var P_factor = 1.0
  let I_ref = 0.2 * cell.param.Q_bol / 3600
  if cell.I_lowpass < -I_ref:
    P_factor = pow(abs(cell.I_lowpass) / I_ref, 1 - cell.param.peukert)
  return cell.param.Q_bol * T_factor * P_factor * cell.soh


proc get_soc(cell: Cell): Soc =
  let Q_effective = cell.get_Q_effective()
  return (Q_effective + cell.C) / Q_effective


proc T_to_R_factor(cell: Cell): float =
  return interpolate(cell.param.T_R_tab, cell.T)


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


proc update_R(cell: Cell) =
  let param = cell.param

  let T_factor = cell.T_to_R_factor()
  let SOH_factor = cell.SOH_to_R_factor()
  let SOC_factor = cell.SOC_to_R_factor()

  cell.model.R0 = param.R0 * T_factor * SOH_factor
  cell.model.R1 = param.R1 * T_factor * SOC_factor * SOH_factor
  cell.model.R2 = param.R2 * T_factor * SOC_factor * SOH_factor
  cell.R = cell.model.R0 + cell.model.R1 + cell.model.R2


proc update_soh(cell: Cell, dt: Interval) =
  let param = cell.param

  # Arrhenius equation to calculate a degradation factor based on the cell's
  # current temperature
  proc arrhenius(ap: ArrheniusParam, T: Temperature): float =
    return ap.A * exp(-ap.Ea / (8.314 * (T + 273.15)))

  cell.dsoh = 0.0

  # 'static' calendar aging
  cell.dsoh += -arrhenius(param.ap_static, cell.T)

  # 'stress' aging due to cycling
  let dC = abs(cell.I) * dt
  if dC > 0:
    # Stress factors: power, SOC window
    let I_nominal = cell.param.Q_bol / 3600
    let power_stress = pow(abs(cell.I) / I_nominal, 1.5)
    let soc_stress = interpolate(param.SOC_stress_Tab, cell.get_soc())
    let degradation_rate = arrhenius(param.ap_stress, cell.T)
    cell.dsoh += -degradation_rate * power_stress * soc_stress
  
  cell.soh += cell.dsoh * dt

  if cell.soh <= 0.0:
    raise newException(ValueError, "Cell SOH dropped to zero")


proc update(cell: Cell, I: Current, dt: Interval) =
  let param = cell.param

  # Update the current in the cell. For the peukert calculation a lowpass
  # filter is used to smooth out short current spikes which can lead to
  # numerical instabilties in parallel cell configurations.
  cell.I = I
  cell.I_lowpass = (cell.I_lowpass * 0.9) + (I * 0.1)

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

  # Update the cell charge, taking into account the charge efficiency
  var dC = I * dt
  if I > 0.0:
    let dynamic_charge_eff = param.charge_eff - (cell.R * param.R_efficiency_factor)
    dC *= dynamic_charge_eff

  # Calculate leak current; given current is at 20°C, adjust for temperature
  let I_leak = param.I_leak_20 * pow(2, (cell.T - 20.0) / 10.0)
  dC += I_leak * dt

  cell.C += dC
  cell.Q_total += abs(dC)

  # Update cell temperature
  let P_R0 = I * I * cell.model.R0
  let P_R1 = I_R1 * I_R1 * cell.model.R1
  let P_R2 = I_R2 * I_R2 * cell.model.R2
  let P_loss = cell.U * max(0, I) * (1.0 - param.charge_eff)
  let P_leak = abs(I_leak) * cell.U
  let P_dis = P_R0 + P_R1 + P_R2 + P_loss + P_leak
  let P_env = (cell.T - cell.T_ambient) / param.Rth
  cell.T += (P_dis - P_env) * dt / param.Cth

  # Update state of health
  cell.update_soh(dt)

  # Get OCV from SOC
  let soc = cell.get_soc()
  let U_ocv = SOC_to_U(param, soc)
  cell.U_src = U_ocv + cell.model.U1 + cell.model.U2
  cell.U = U_ocv + cell.model.U0 + cell.model.U1 + cell.model.U2


proc newCell(sim: Simulation, param: CellParam): Cell =
  var cell = new Cell
  cell.id = sim.cells.len()
  cell.sim = sim
  cell.param = param
  cell.model = CellModel()
  cell.C = 0.0
  cell.R = cell.param.R0 + cell.param.R1 + cell.param.R2
  cell.T = 20.0
  cell.T_ambient = 20.0
  cell.soh = 1.0

  let fname = fmt"/tmp/cell_{cell.id:02}.log"
  cell.fd_log = open(fname, fmWrite)

  # Run one step to initialize cell state
  cell.update_R()
  cell.update(0.0, 0.0)

  sim.cells.add(cell)

  # Deviations
  cell.param.Q_bol *= rand(0.9 .. 1.00)
  cell.param.R0 *= rand(0.95 .. 1.00)

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

  proc gen_graph(col: int, ylabel: string, pres: openArray[string]) =
    var ts: seq[string]
    for cell in sim.cells.mitems:
      ts.add(&""" "/tmp/cell_{cell.id:02}.log" u {col} w l""")

    l("")
    for pre in pres:
      l(pre)
    l(&"""set ylabel "{ylabel}"""")
    l(&"""plot {ts.join(", ")}""")

  gen_graph(1, "I (A)",    [ "unset yrange" ])
  gen_graph(2, "U (V)",    [ "set yrange [2.3:4.4]" ])
  gen_graph(3, "SOC (%)",  [ "set yrange [-0.1:1.1]" ])
  gen_graph(4, "T (°C)",   [ "unset yrange" ])
  gen_graph(5, "SOH (%)",  [ "set yrange [-0.1:1.1]" ])
  #gen_graph(6, "dSOH",     [ "unset yrange" ])

  fd.write("unset multiplot\n")
  fd.close()


proc report(cell: Cell) =
  if cell.sim.cycle_number mod cell.sim.report_every_n == 0:
    let line = fmt"{cell.I:>4.2f} {cell.U:>6.3f} {cell.get_soc():>4.3f} {cell.T:>4.2f} {cell.soh:>4.2f} {log10(-cell.dsoh):>g}"
    cell.fd_log.writeLine(line)


proc balance(sim: Simulation, pack: var Pack) =
  var U_min = 1e6
  for module in pack.modules:
    for cell in module.cells:
      if cell.U < U_min:
        U_min = cell.U
  for module in pack.modules.mitems:
    var I_balance = 0.0
    for cell in module.cells:
      let dU = cell.U - U_min
      if dU > 0.01:
        I_balance -= pack.balancer.I
    module.I_balance = I_balance


proc step(sim: Simulation, I_pack: Current): Voltage =

  sim.balance(sim.pack)

  var i = 0

  for module in sim.pack.modules:
   
    # Current in the module is pack current + balancing current
    let I_module = I_pack + module.I_balance

    var sum_U_div_R = 0.0
    var sum_1_div_R = 0.0

    # Calculate parallel resistance of all cells in the module
    for cell in module.cells:
      cell.update_R()
      sum_U_div_R += cell.U_src / cell.model.R0
      sum_1_div_R += 1.0 / cell.model.R0

    let U_module = (I_module + sum_U_div_R) / sum_1_div_R
    result += U_module

    for cell in module.cells:
      cell.report()
      inc i
      let I_cell = (U_module - cell.U_src) / cell.model.R0
      cell.update(I_cell, sim.dt)
        
  sim.time += sim.dt


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
    var U = sim.step(I)
    if U < U_min:
      return
    for cell in sim.cells:
      if cell.U < 2.5:
        return

proc charge(sim: Simulation, I: Current, U_max: Voltage) =
  while true:
    var U_pack = sim.step(I)
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
    var U_pack = sim.step(clamp(I_pack, 0.0, I_set))

    let err = U_set - U_pack
    err_int += err * sim.dt
    err_int = clamp(err_int, -2.0, 2.0)
  
    I_pack = (kP * err) + (kI * err_int)

    if sim.time > t_max:
      raise newException(ValueError, "CC/CV charge timeout")
    

proc sleep(sim: Simulation, d: Duration) =
  let t_end = sim.time + d
  while sim.time < t_end:
    discard sim.step(0.0)



let param = CellParam(
  vendor: "Samsung",
  model:  "INR18650-32E",
  R0:     0.025,
  R1:     0.015,
  R2:     0.010,
  C1:  4000.00,
  C2: 30000.00,
  Cth:   25.00,
  Rth:    5.00,
  Q_bol: Q_from_Ah(3.2),
  I_leak_20: -1.4e-3,
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
  charge_eff: 0.96,
  peukert: 1.03,
  R_efficiency_factor: 5.0,
  ap_static: ArrheniusParam(
    A: 200.0,
    Ea: 55.0e3,
  ),
  ap_stress: ArrheniusParam(
    A: 300,
    Ea: 54.0e3,
  )
)



proc test1(sim: Simulation) =
  sim.sleep(600)
  sim.discharge(-8.0, sim.pack.U_empty)
  sim.sleep(1200)
  #sim.charge(+4.0, sim.pack.U_full)
  sim.charge_CC_CV(+4.0, sim.pack.U_full)
  sim.sleep(600)

proc test2(sim: Simulation) =
  sim.sleep(4 * 3600)


proc test3(sim: Simulation) =
  proc drive(sim: Simulation, I: Current, d: Duration) =
    let t_end = sim.time + d
    while sim.time < t_end:
      discard sim.step(I)
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


proc run(sim: Simulation, fn: proc(sim: Simulation), count: int, n_report: int=5) =
  sim.report_every_n = max(1, count div n_report)
  for i in 0 ..< count:
    sim.cycle_number = i
    fn(sim)
  let t = sim.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  echo fmt"Completed {count} cycles, total time {days}d {hours}h {minutes}m"
  echo fmt"SOH range: {sim.cells.mapIt(it.soh).min:.3f} .. {sim.cells.mapIt(it.soh).max:.3f}"


var sim = newSimulation(5.0)
sim.pack = sim.newPack(n_series=2, n_parallel=2, param)
sim.pack.balancer.I = 0.200

sim.cells[0].param.R0 *= 1.2
sim.cells[0].param.Q_bol *= 1.0

sim.gen_gnuplot("view.gp")
sim.run(test1, count=4, n_report=10)


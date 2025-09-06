
import std / [math, strformat, random, os]

import types
import rc
import cell
import pack
import battery
import model
import misc
import tests_eis
import tests_commute
import gnuplot

let param = CellParam(
  vendor:   "Samsung",
  model:    "INR18650-32E",
  RC_dc:    RCParam(R: 0.040),
  RC_trans: RCParam(R: 0.015, C:  4000.0),
  RC_diff: @[ 
    RCParam(R: 0.008, C:    11_200),
    RCParam(R: 0.006, C:    50_000),
    RCParam(R: 0.004, C:   240_000),
    RCParam(R: 0.002, C: 1_200_000),
  ],
  RCt_core: RCtParam(R: 2.5, C:   150.0),
  RCt_cell: RCtParam(R: 5.0, C:    30.0),
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
  # Samsung-INR18650-32E.pdf 7.5
  T_cap_tab: @[
    (-10.0, 0.60),
    ( 25.0, 1.00),
    ( 40.0, 1.00)
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
  # energies-12-02685.pdf, figure 5
  entropy_tab: @[
    (0.0,  0.0002),
    (0.2,  0.0001),
    (0.5, -0.0001),
    (0.8, -0.0003),
    (1.0,  0.0004)
  ],
  # Samsung-INR18650-32E.pdf 7.5
  charge_eff: 0.97,
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

let batt_param = BatteryParam(
  RCt_air: RCtParam(R: 2.0, C: 10.0),
  RCt_case: RCtParam(R: 1.0, C: 250.0),
  T_env: 20.0,
)

proc test_cycle(sim: Simulation) =
  sim.sleep(600)
  sim.discharge(-5.6, sim.battery.pack.U_empty)
  sim.sleep(3600)
  #sim.charge(+4.0, sim.battery.pack.U_full)
  sim.charge_CC_CV(+4.0, sim.battery.pack.U_full)
  sim.sleep(600)


proc test_sleep(sim: Simulation) =
  sim.sleep(24 * 3600)




var sim = newSimulation(5.0)
sim.battery.init(batt_param)
sim.battery.pack.init(n_series=4, n_parallel=4, param)
sim.battery.balancer.I = 0.200

sim.run(test_cycle, count=1, n_report=1)
#sim.run(test_EIS)
#sim.run(test_commute)
sim.gen_gnuplot("battery.gp")


import std / [math]

import types
import rc
import cell
import pack
import battery
import balancer
import model
import misc
import tests_eis
import tests_commute
import gnuplot

let cell_param = CellParam(
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
  SOC_to_U: mkLut @[
    (0.00, 2.50),
    (0.05, 2.91),
    (0.10, 3.12),
    (0.15, 3.25),
    (0.20, 3.34),
    (0.25, 3.42),
    (0.30, 3.50),
    (0.35, 3.57),
    (0.40, 3.62),
    (0.45, 3.68),
    (0.50, 3.72),
    (0.55, 3.77),
    (0.60, 3.81),
    (0.65, 3.86),
    (0.70, 3.92),
    (0.75, 3.98),
    (0.80, 4.03),
    (0.85, 4.07),
    (0.90, 4.09),
    (0.95, 4.12),
    (1.00, 4.20),
  ],
  # Samsung-INR18650-32E.pdf 7.5
  T_to_cap: mkLut @[
    (-10.0, 0.60),
    ( 25.0, 1.00),
    ( 40.0, 1.00)
  ],
  T_to_R: mkLut @[
    (-20.0, 3.0),
    (  0.0, 1.8),
    ( 25.0, 1.0)
  ],
  SOH_to_R: mkLut @[
    (0.0, 1.5),
    (0.2, 1.3),
    (0.4, 1.2),
    (0.6, 1.1),
    (0.8, 1.05),
    (1.0, 1.0)
  ],
  # EVE INR21700-50E Test Report
  SOC_to_R: mkLut @[
    (0.0, 3.95),
    (0.1, 1.52),
    (0.2, 1.11),
    (0.3, 1.02),
    (0.4, 1.00),
    (0.5, 1.01),
    (0.6, 1.06),
    (0.7, 1.05),
    (0.8, 1.05),
    (0.9, 1.04),
    (1.0, 1.29),
  ],
  SOC_to_stress: mkLut @[
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
  SOC_to_entropy: mkLut @[
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

let balancer_param = BalancerParam(
  I: -0.100,
  U_min: 3.60,
  U_max: 4.30,
  U_delta: 0.02
)

let batt_param = BatteryParam(
  n_series: 4,
  n_parallel: 4,
  RCt_air: RCtParam(R: 2.0, C: 10.0),
  RCt_case: RCtParam(R: 1.0, C: 250.0),
  T_env: 20.0,
  cell_param: cell_param,
  balancer_param: balancer_param,
)


proc test_cycle(model: Model) =
  model.sleep(1800)
  model.discharge(-5.6, model.battery.pack.U_empty)
  model.sleep(3600)
  #model.charge(+4.0, model.battery.pack.U_full)
  model.charge_CC_CV(+4.0, model.battery.pack.U_full)
  model.sleep(1800)


proc test_sleep(model: Model) =
  model.sleep(24 * 3600)


block:
  var model = newModel(dt=5.0)
  model.battery.init(batt_param)

  model.run(test_cycle, count=1, n_report=10)
  #model.run(test_EIS)
  #model.run(test_commute)
  model.gen_gnuplot("battery.gp")


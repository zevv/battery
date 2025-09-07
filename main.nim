
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
    (0.0, 2.500),
    (0.02, 2.710),
    (0.04, 2.868),
    (0.06, 2.972),
    (0.08, 3.053),
    (0.1, 3.115),
    (0.12, 3.168),
    (0.14, 3.212),
    (0.16, 3.258),
    (0.18, 3.304),
    (0.2, 3.347),
    (0.22, 3.386),
    (0.24, 3.422),
    (0.26, 3.460),
    (0.28, 3.484),
    (0.3, 3.500),
    (0.32, 3.519),
    (0.34, 3.545),
    (0.36, 3.571),
    (0.38, 3.594),
    (0.4, 3.615),
    (0.42, 3.638),
    (0.44, 3.659),
    (0.46, 3.679),
    (0.48, 3.700),
    (0.5, 3.722),
    (0.52, 3.744),
    (0.54, 3.766),
    (0.56, 3.786),
    (0.58, 3.805),
    (0.6, 3.823),
    (0.62, 3.839),
    (0.64, 3.855),
    (0.66, 3.873),
    (0.68, 3.893),
    (0.7, 3.913),
    (0.72, 3.932),
    (0.74, 3.953),
    (0.76, 3.977),
    (0.78, 4.005),
    (0.8, 4.032),
    (0.82, 4.055),
    (0.84, 4.072),
    (0.86, 4.081),
    (0.88, 4.086),
    (0.9, 4.090),
    (0.92, 4.094),
    (0.94, 4.100),
    (0.96, 4.120),
    (0.98, 4.150),
    (1.0, 4.200),
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


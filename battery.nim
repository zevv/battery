
import types
import rc
import cell
import pack
import balancer

type 

  BatteryParam* = object
    RCt_air*: RCtParam # air to case
    RCt_case*: RCtParam # case to environment
    T_env*: Temperature # ambient temperature
    n_series*: int
    n_parallel*: int
    cell_param*: CellParam
    balancer_param*: BalancerParam

  Battery* = object
    pack*: Pack
    param*: BatteryParam
    RCt_air*: RCtModel
    RCt_case*: RCtModel
    balancer*: Balancer
    dis_FET_on*: bool
    chg_FET_on*: bool
    # outputs
    U*: Voltage
    U_cell_min*: Voltage
    U_cell_max*: Voltage


proc update_temperature(battery: var Battery, dt: Interval) =
  var P_pack = battery.pack.get_P_heat()

  battery.RCt_air.update(battery.param.RCt_air, 
             P_pack, battery.RCt_case.T, dt)
  battery.RCt_case.update(battery.param.RCt_case, 
             battery.RCt_air.P, battery.param.T_env, dt)


proc step*(battery: var Battery, I: Current, dt: Interval) =

  var I = I
  if I > 0 and not battery.chg_FET_on:
    I = 0.0
  if I < 0 and not battery.dis_FET_on:
    I = 0.0
 
  let U_pack = battery.pack.step(I, battery.RCt_air.T, dt)

  battery.update_temperature(dt)

  let Us = battery.pack.get_U_cells()
  let Is = battery.balancer.step(I, Us)
  battery.pack.set_I_balance(Is)

  battery.U = U_pack
  battery.U_cell_min = Us.min()
  battery.U_cell_max = Us.max()

  if battery.dis_FET_on and battery.U_cell_min < 2.50:
    echo "* Switching off discharge FET, U_cell_min = ", battery.U_cell_min
    battery.dis_FET_on = false

  if battery.chg_FET_on and battery.U_cell_max > 4.25:
    echo "* Switching off charge FET, U_cell_max = ", battery.U_cell_max
    battery.chg_FET_on = false


proc reset*(battery: var Battery) =
  battery.dis_FET_on = true
  battery.chg_FET_on = true


proc report*(battery: Battery, t: Interval) =
  battery.pack.report(t, battery.RCt_case.T)


proc init*(battery: var Battery, param: BatteryParam) =
  battery.param = param
  battery.RCt_case.T = 20.0
  battery.RCt_air.T = 20.0

  battery.dis_FET_on = true
  battery.chg_FET_on = true

  battery.pack.init(param.n_series, param.n_parallel, param.cell_param)
  battery.balancer.init(param.n_series, param.balancer_param)

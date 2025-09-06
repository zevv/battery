
import types
import rc
import pack
import balancer

type 

  BatteryParam* = object
    RCt_case*: RCtParam # case to environment
    RCt_air*: RCtParam # case to environment
    T_env*: Temperature # ambient temperature

  Battery* = object
    pack*: Pack
    param*: BatteryParam
    RCt_case*: RCtModel
    RCt_air*: RCtModel
    balancer*: Balancer


proc update_temperature(battery: var Battery, dt: Interval) =
  var P_pack = battery.pack.get_P_heat()

  battery.RCt_air.update(battery.param.RCt_air, 
             P_pack, battery.RCt_case.T, dt)
  battery.RCt_case.update(battery.param.RCt_case, 
             battery.RCt_air.P, battery.param.T_env, dt)


proc step*(battery: var Battery, I: Current, dt: Interval): Voltage =
 
  let U_pack = battery.pack.step(I, battery.RCt_air.T, dt)

  battery.update_temperature(dt)

  let Us = battery.pack.get_U_cells()
  let Is = battery.balancer.step(Us)
  battery.pack.set_I_balance(Is)

  U_pack


proc report*(battery: Battery, t: Interval) =
  battery.pack.report(t, battery.RCt_case.T)


proc init*(battery: var Battery, param: BatteryParam) =
  battery.param = param
  battery.RCt_case.T = 20.0

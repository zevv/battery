
import types
import rc
import cell
import pack


type 

  Balancer = object
    I*: Current

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
  var P_cells = 0.0
  for module in battery.pack.modules.mitems:
    for cell in module.cells.mitems:
      P_cells += cell.RCt_cell.P

  battery.RCt_air.update(battery.param.RCt_air, 
             P_cells, battery.RCt_case.T, dt)
  battery.RCt_case.update(battery.param.RCt_case, 
             battery.RCt_air.P, battery.param.T_env, dt)


proc step*(battery: var Battery, I: Current, dt: Interval): Voltage =
  let U_pack = battery.pack.step(I, battery.RCt_air.T, dt)
  battery.update_temperature(dt)
  U_pack


proc report*(battery: Battery, t: Interval) =
  battery.pack.report(t, battery.RCt_case.T)


proc init*(battery: var Battery, param: BatteryParam) =
  battery.param = param
  battery.RCt_case.T = 20.0


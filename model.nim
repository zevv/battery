
import std/[sequtils, strformat, math]

import types
import pack
import battery

type 

  Model* = ref object
    battery*: Battery
    time*: float
    time_report: float
    steps: int
    dt*: float
    cycle_number: int
    report_every_n: int


proc step*(model: Model, I: Current, dt: Interval): Voltage =

  let U_batt = model.battery.step(I, dt)

  if model.cycle_number mod model.report_every_n == 0:
    model.battery.report(model.time_report)
    model.time_report += dt
  
  model.time += dt
  inc model.steps

  U_batt


proc newModel*(dt: Interval): Model =
  result = new Model
  result.dt = dt


proc run*(model: Model, fn: proc(model: Model), count: int=1, n_report: int=1) =

  model.report_every_n = max(1, count div n_report)

  for i in 0 ..< count:
    model.cycle_number = i
    fn(model)
  let t = model.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  stderr.write &"Completed {count} cycles, {model.steps} steps, total time {days}d {hours}h {minutes}m\n"



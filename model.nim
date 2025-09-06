
import std/[sequtils, strformat, math, times]

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

  model.battery.step(I, dt)

  if model.cycle_number mod model.report_every_n == 0:
    model.battery.report(model.time_report)
    model.time_report += dt
  
  model.time += dt
  inc model.steps

  model.battery.U


proc newModel*(dt: Interval): Model =
  result = new Model
  result.dt = dt


proc run*(model: Model, fn: proc(model: Model), count: int=1, n_report: int=1) =

  model.report_every_n = max(1, count div n_report)

  let t_start = epochTime()
  for i in 0 ..< count:
    model.cycle_number = i
    model.battery.reset()
    fn(model)
  let t_stop = epochTime()
  let duration = t_stop - t_start

  let t = model.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  stderr.write &"* Completed {count} cycles in {model.steps} steps, simulation time {days}d {hours}h {minutes}m\n"
  stderr.write &"* Real time: {duration:.2f} seconds, speedup {model.time / duration:.3g}x, steps/s: {model.steps.float / duration / 1000.0:.1f}k\n"




import std/[sequtils, tables, strutils, times, os, strformat]
import std/[math, random, algorithm, tempfiles, complex]

import types
import rc
import cell
import pack
import battery

type 

  Simulation* = ref object
    pre_hook: proc(sim: Simulation, I_pack: Current)
    battery*: Battery
    time*: float
    time_report: float
    steps: int
    dt*: float
    cycle_number: int
    report_every_n: int


proc step*(sim: Simulation, I: Current, dt: Interval): Voltage =

  if sim.pre_hook != nil:
    sim.pre_hook(sim, I)

  let U_batt = sim.battery.step(I, dt)

  if sim.cycle_number mod sim.report_every_n == 0:
    sim.battery.report(sim.time_report)
    sim.time_report += dt

  sim.time += dt
  inc sim.steps

  U_batt


proc newSimulation*(dt: Interval): Simulation =
  result = new Simulation
  result.dt = dt


proc balance(sim: Simulation, I_pack: Current) =
  if I_pack > 0.0:
  
    let U_min = sim.battery.pack.modules.mapIt(it.U).min
    for module in sim.battery.pack.modules.mitems:
      let dU = module.U - U_min
      if module.U > 4.1 and dU > 0.02:
        module.I_balance = - sim.battery.balancer.I
      else:
        module.I_balance = 0.0


proc run*(sim: Simulation, fn: proc(sim: Simulation), count: int=1, n_report: int=1) =

  sim.pre_hook = balance
  sim.report_every_n = max(1, count div n_report)

  for i in 0 ..< count:
    sim.cycle_number = i
    fn(sim)
  let t = sim.time.int
  let days = t div 86400
  let hours = (t mod 86400) div 3600
  let minutes = (t mod 3600) div 60
  stderr.write &"Completed {count} cycles, {sim.steps} steps, total time {days}d {hours}h {minutes}m\n"



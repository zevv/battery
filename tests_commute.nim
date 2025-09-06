
import std/[ random, math ]

import types
import model
import misc

proc test_commute*(sim: Simulation) =
  proc drive(sim: Simulation, I: Current, d: Duration) =
    let t_end = sim.time + d
    while sim.time < t_end:
      discard sim.step(I, sim.dt)
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
  sim.charge_CC_CV(+4.0, sim.battery.pack.U_full)
  echo "--- Overnight Rest ---"
  let total_cycle_duration = 24.0 * 3600.0
  let time_into_this_cycle = sim.time mod total_cycle_duration
  let remaining_time = total_cycle_duration - time_into_this_cycle
  if remaining_time > 0:
    sim.sleep(remaining_time)


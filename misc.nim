
import model

proc discharge*(sim: Simulation, I: Current, U_min: Voltage) =
  while true:
    var U = sim.step(I, sim.dt)
    if U < U_min:
      return
    for module in sim.pack.modules:
      for cell in module.cells:
        if cell.U < 2.5:
          return

proc charge*(sim: Simulation, I: Current, U_max: Voltage) =
  while true:
    var U_pack = sim.step(I, sim.dt)
    if U_pack > U_max:
      break


proc charge_CC_CV*(sim: Simulation, I_set: Current, U_set: Voltage) =
  
  # PID controller constants for voltage regulation
  let kP = 0.3
  let kI = 2.5

  let t_max = sim.time + 6 * 3600
  var I_pack = I_set
  var err_int = 0.0

  while I_pack > I_set * 0.05:
    let I = clamp(I_pack, 0.0, I_set)
    var U_pack = sim.step(I, sim.dt)

    let err = U_set - U_pack
    err_int += err * sim.dt
    err_int = clamp(err_int, -2.0, 2.0)
  
    I_pack = (kP * err) + (kI * err_int)

    if sim.time > t_max:
      raise newException(ValueError, "CC/CV charge timeout")
    

proc sleep*(sim: Simulation, d: Duration) =
  let t_end = sim.time + d
  while sim.time < t_end:
    discard sim.step(0.0, sim.dt)





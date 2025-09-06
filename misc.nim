
import types
import model

proc discharge*(model: Model, I: Current, U_min: Voltage) =
  while true:
    var U = model.step(I, model.dt)
    if U < U_min:
      return
    for module in model.battery.pack.modules:
      for cell in module.cells:
        if cell.U < 2.5:
          return

proc charge*(model: Model, I: Current, U_max: Voltage) =
  while true:
    var U_pack = model.step(I, model.dt)
    if U_pack > U_max:
      break


proc charge_CC_CV*(model: Model, I_set: Current, U_set: Voltage) =
  
  # PID controller constants for voltage regulation
  let kP = 0.3
  let kI = 2.5

  let t_max = model.time + 6 * 3600
  var I_pack = I_set
  var err_int = 0.0

  while I_pack > I_set * 0.05:
    let I = clamp(I_pack, 0.0, I_set)
    var U_pack = model.step(I, model.dt)

    let err = U_set - U_pack
    err_int += err * model.dt
    err_int = clamp(err_int, -2.0, 2.0)
  
    I_pack = (kP * err) + (kI * err_int)

    if model.time > t_max:
      raise newException(ValueError, "CC/CV charge timeout")
    

proc sleep*(model: Model, d: Duration) =
  let t_end = model.time + d
  while model.time < t_end:
    discard model.step(0.0, model.dt)





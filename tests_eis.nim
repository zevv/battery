
import std / [math, strformat]

import model

proc test_EIS_f(sim: Simulation, freq: float, I_bias: Current) =
  let cycles = 3
  let steps_per_cycle = 100
  let steps = cycles * steps_per_cycle
  let dt = 1.0 / (freq * steps_per_cycle.float)
  stderr.write(&"EIS {freq:>8.3f} Hz, dt={dt:>6.4g} s, steps={steps}\n")
  let I_amp = +0.020
  var Zr = 0.0
  var Zi = 0.0
  var n = 0
  let t_start = sim.time
  while n < steps:
    let t = sim.time - t_start
    let ref_sin = sin(TAU * freq * t)
    let ref_cos = cos(TAU * freq * t)
    var U = sim.step(ref_sin * I_amp + I_bias, dt)
    Zr += U * ref_sin
    Zi += U * ref_cos
    inc n
  let scale = 2.0 / (I_amp * steps.float)
  Zr *= scale
  Zi *= scale
  let pha = arctan2(Zi, Zr) * 180.0 / PI
  let mag = sqrt(Zr*Zr + Zi*Zi)
  echo freq, " ", Zr, " ", Zi, " ", mag, " ", pha


#proc test_EIS_f2(sim: Simulation, freq: float, I_bias: Current) =
#  let w = TAU * freq
#  var Z = complex(param.RC_dc.R)
#  Z += param.RC_trans.R / complex(1.0, w * param.RC_trans.R * param.RC_trans.C)
#  for rc in param.RC_diff:
#    Z += rc.R / complex(1.0, w * rc.R * rc.C)
#  let mag = abs(Z)
#  let pha = arctan2(Z.im, Z.re) * 180.0 / PI
#  echo freq, " ", Z.re, " ", Z.im, " ", mag, " ", pha
#

# https://pubs.acs.org/doi/pdf/10.1021/acsmeasuresciau.2c00070?ref=article_openPDF

proc discharge_time(sim: Simulation, I: Current, t: Duration) =
  let t_start = sim.time
  while sim.time < t_start + t:
    discard sim.step(I, sim.dt)


proc test_EIS*(sim: Simulation) =
  let I_bias = -0.001
  var f = 0.0003
  sim.discharge_time(I_bias, 180)
  while f < 100:
    #sim.charge_CC_CV(+4.0, sim.pack.U_full)
    sim.test_EIS_f(f, I_bias)
    f *= 1.5



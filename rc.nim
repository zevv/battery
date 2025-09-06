
import types

type 
  
  RCParam* = object
    R*: Resistance # Ω for electrical
    C*: Capacitance # F for electrical

  RCtParam* = object
    R*: Resistance # K/W for thermal
    C*: Capacitance # J/K for thermal
  
  RCModel* = object
    R*: Resistance
    U*: Voltage
    I_R*: Current
    I_C: Current

  RCtModel* = object
    T*: Temperature
    P*: Power



proc update*(m: var RCtModel, rp: RCtParam, P_in: Power, T_out: Temperature, dt: Interval) =
  m.P = (m.T - T_out) / rp.R
  m.T += (P_in - m.P) * dt / rp.C
  assert(m.T > -20)
  assert(m.T <  80)


# Timestep the RC equivalent circuit model. Special case for C=0 (pure resistor).

proc update*(rc: var RCModel, rp: RCParam, I: Current, dt: Interval) =
  if rp.C > 0.0:
    rc.I_R = rc.U / rp.R
    rc.I_C = I - rc.I_R
    rc.U += dt * rc.I_C / rp.C
  else:
    rc.I_R = I
    rc.I_C = 0.0
    rc.U = I * rc.R



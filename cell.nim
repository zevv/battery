
import std/[sequtils, tables, strutils, times, os, strformat]
import std/[math, random, algorithm, tempfiles, complex]

import types
import rc

type
  
  ArrheniusParam* = object
    A*: float # pre-exponential factor, 1/s
    Ea*: float # activation energy, J/mol

  CellParam* = object
    vendor*: string
    model*: string
    RC_dc*: RCParam # DC resistance
    RC_trans*: RCParam # charge transfer
    RC_diff*: seq[RCParam] # diffusion model
    Q_bol*: Charge # nominal capacity at 1C, 3600*Ah
    I_leak_20*: Current # self-discharge current at 20°C
    soc_tab*: SocTab # OCV vs SOC
    T_cap_tab*: Tab[Temperature, float] # capacity factor vs temperature
    T_R_tab*: Tab[Temperature, float] # resistance factor vs temperature
    SOH_R_tab*: Tab[Soh, float] # resistance factor vs SOH
    SOC_R_tab*: Tab[Soc, float] # resistance factor vs SOC
    SOC_stress_Tab*: Tab[Soc, float] # 'stress' factor vs SOC
    entropy_tab*: Tab[Soc, float] # entropy coefficient vs SOC
    RCt_core*: RCtParam # RC thermal model core to case
    RCt_cell*: RCtParam # RC thermal model case to environment
    charge_eff*: float # nominal charge efficiency
    peukert*: float # Peukert exponent
    R_efficiency_factor*: float # approximates charge efficiency drop for various chemical effects
    ap_static*: ArrheniusParam # calendar aging
    ap_stress*: ArrheniusParam # cycling aging

  Cell* = object
    param*: CellParam
    fd_log*: File
    RC_dc*: RCModel # 
    RC_trans*: RCModel # charge transfer
    RC_diff*: seq[RCModel] # diffusion
    Q*: Charge # energy taken out
    R*: Resistance
    T_env*: Temperature
    RCt_core*: RCtModel
    RCt_cell*: RCtModel
    I_leak*: Current
    I*: Current
    I_lowpass*: Current
    soc*: Soc
    U*: Voltage # terminal voltage
    U_src*: Voltage # source voltage (without R0 drop)
    soh*: Soh


proc SOC_to_U(cp: CellParam, soc: Soc): float =
  let n = len(cp.soc_tab)
  if soc <= 0.0:
    return cp.soc_tab[0]
  if soc >= 1.0:
    return cp.soc_tab[n - 1]
  let idx = int(soc * (n - 1).float)
  let f1 = cp.soc_tab[idx]
  let f2 = cp.soc_tab[idx + 1]
  let soc1 = float(idx) / (n - 1).float
  let soc2 = float(idx + 1) / (n - 1).float
  return f1 + (f2 - f1) * (soc - soc1) / (soc2 - soc1)


proc update_soc(cell: var Cell) =

  # Temperature factor
  let T_factor = interpolate(cell.param.T_cap_tab, cell.RCt_core.T)

  # Peukert factor
  var P_factor = 1.0
  let I_ref = 0.2 * cell.param.Q_bol / 3600
  if cell.I_lowpass < -I_ref:
    P_factor = pow(abs(cell.I_lowpass) / I_ref, 1 - cell.param.peukert)
  let Q_effective = cell.param.Q_bol * T_factor * P_factor * cell.soh

  cell.soc = (Q_effective + cell.Q) / Q_effective


proc update_R*(cell: var Cell) =
  let param = cell.param

  let T_factor = interpolate(param.T_R_tab, cell.RCt_core.T)
  let SOH_factor = interpolate(param.SOH_R_tab, cell.soh)
  let SOC_factor = interpolate(param.SOC_R_tab, cell.soc)

  cell.RC_dc.R = param.RC_dc.R * T_factor * SOH_factor
  cell.RC_trans.R = param.RC_trans.R * T_factor * SOC_factor * SOH_factor

  for i, rc in cell.RC_diff.mpairs:
    rc.R = param.RC_diff[i].R * T_factor * SOC_factor * SOH_factor


# https://en.wikipedia.org/wiki/Arrhenius_equation
proc arrhenius(ap: ArrheniusParam, T: Temperature): float =
  return ap.A * exp(-ap.Ea / (8.314 * (T + 273.15)))


proc update_soh(cell: var Cell, dt: Interval) =
  let param = cell.param

  var soh_rate = 0.0

  # 'static' calendar aging
  soh_rate += -arrhenius(param.ap_static, cell.RCt_core.T)

  # Stress factors: power, SOC window
  let I_nominal = cell.param.Q_bol / 3600
  let power_stress = pow(abs(cell.I) / I_nominal, 1.5)
  let soc_stress = interpolate(param.SOC_stress_Tab, cell.soc)
  let degradation_rate = arrhenius(param.ap_stress, cell.RCt_core.T)
  soh_rate += -degradation_rate * power_stress * soc_stress
 
  cell.soh += soh_rate * dt

  if cell.soh <= 0.0:
    raise newException(ValueError, "Cell SOH dropped to zero")


proc update_temperature(cell: var Cell, I: Current, T_env: Temperature, dt: Interval) =

  let param = cell.param

  # Power losses
  let P_R0 = I * I * cell.RC_dc.R
  let P_R_trans = cell.RC_trans.I_R * cell.RC_trans.I_R * cell.RC_trans.R
  let P_R_diff = cell.RC_diff.foldl(a + b.I_R * b.I_R * b.R, 0.0)
  let P_loss = cell.U * max(0, I) * (1.0 - param.charge_eff)
  let P_leak = abs(cell.I_leak) * cell.U
  let P_rev = interpolate(param.entropy_tab, cell.soc) * (cell.RCt_core.T+273.15) * I
  let P_dis = P_R0 + P_R_trans + P_R_diff + P_loss + P_leak + P_rev

  cell.RCt_core.update(param.RCt_core, P_dis, cell.RCt_cell.T, dt)
  cell.RCt_cell.update(param.RCt_cell, cell.RCt_core.P, T_env, dt)


proc update_charge(cell: var Cell, I: Current, dt: Interval) =
  let param = cell.param

  # Update cell charge, taking into account the charge efficiency
  var dQ = I * dt
  if I > 0.0:
    let dynamic_charge_eff = param.charge_eff - (cell.RC_dc.R * param.R_efficiency_factor)
    dQ *= dynamic_charge_eff

  # Calculate leak current; given current is at 20°C, adjust for temperature
  cell.I_leak = param.I_leak_20 * pow(2, (cell.RCt_core.T - 20.0) / 10.0)
  dQ += cell.I_leak * dt

  cell.Q += dQ


proc update_voltage(cell: var Cell) =
  let param = cell.param
  let soc = cell.soc
  let U_ocv = SOC_to_U(param, soc)
  let U_diff = cell.RC_diff.foldl(a + b.U, 0.0)
  cell.U_src = U_ocv + cell.RC_trans.U + U_diff
  cell.U = U_ocv + cell.RC_dc.U + cell.RC_trans.U + U_diff


proc update*(cell: var Cell, I: Current, T_env: Temperature, dt: Interval) =
  let param = cell.param

  cell.I = I
  cell.I_lowpass = (cell.I_lowpass * 0.9) + (I * 0.1)

  cell.RC_dc.update(param.RC_dc, I, dt)
  cell.RC_trans.update(param.RC_trans, I, dt)
  for i, rc in cell.RC_diff.mpairs:
    rc.update(param.RC_diff[i], I, dt)

  cell.update_soc()
  cell.update_voltage()
  cell.update_charge(I, dt)
  cell.update_temperature(I, T_env, dt)
  cell.update_soh(dt)


proc report*(cell: var Cell, time: float, T_case: Temperature) =
  let line = &"{time / 60} {cell.I:>4.2f} {cell.U:>6.3f} {cell.soc:>4.3f} {cell.RCt_core.T:>5.3f} {T_case:>5.3f} {cell.soh:>4.2f}"
  cell.fd_log.writeLine(line)

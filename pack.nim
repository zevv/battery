
import types
import cell

type 
  
  Module* = object
    I_balance*: Current
    U*: Voltage
    cells*: seq[Cell]
  
  Pack* = object
    I*: Current
    U_empty*: Voltage
    U_full*: Voltage
    modules*: seq[Module]


proc init*(pack: var Pack, n_series: int, n_parallel: int, param: CellParam) =

  pack.modules = newSeq[Module](n_series)
  for module in pack.modules.mitems:
    module.cells = newSeq[Cell](n_parallel)
    for cell in module.cells.mitems:
      cell.init(param)

  pack.U_empty = n_series.float * 2.50
  pack.U_full = n_series.float * 4.20


proc step*(pack: var Pack, I_pack: Current, T: Temperature, dt: Interval): Voltage =

  for module in pack.modules.mitems:
   
    # Calculate parallel resistance of all cells in the module
    var sum_U_div_R = 0.0
    var sum_1_div_R = 0.0
    for cell in module.cells.mitems:
      cell.update_R()
      sum_U_div_R += cell.U_src / cell.RC_dc.R
      sum_1_div_R += 1.0 / cell.RC_dc.R

    # Current in the module is pack current + balancing current
    let I_module = I_pack + module.I_balance
    module.U = (I_module + sum_U_div_R) / sum_1_div_R
    result += module.U

    # Update each cell in the module with the calculated cell current
    for cell in module.cells.mitems:
      let I_cell = (module.U - cell.U_src) / cell.RC_dc.R
      cell.step(I_cell, T, dt)


proc report*(pack: Pack, t: Interval, T_env: Temperature) =
  for module in pack.modules:
    for cell in module.cells:
      cell.report(t, T_env)


import types

type

  BalancerParam* = object
    I*: Current # balancing current
    U_min*: Voltage # voltage to start balancing
    U_max*: Voltage # voltage to alays balance
    U_delta*: Voltage

  Balancer* = object
    param*: BalancerParam
    I_module*: seq[Current] # balancing currents for each cell


proc init*(balancer: var Balancer, n_modules: int, param: BalancerParam) =
  balancer.param = param
  balancer.I_module = newSeq[Current](n_modules)


proc step*(balancer: var Balancer, I: Current, U_cells: seq[Voltage]): seq[Current] =
  let param = balancer.param

  if I > 0.0:
    let U_min = U_cells.min()

    for i, U in U_cells:
      let dU = U - U_min
      var I = 0.0
      if (U > param.U_min and dU > param.U_delta) or U > param.U_max:
        I = param.I
      balancer.I_module[i] = I

  balancer.I_module

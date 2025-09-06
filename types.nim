

type 

  Soc* = float
  Soh* = float
  Temperature* = float
  Resistance* = float
  Capacitance* = float
  Voltage* = float
  Current* = float
  Charge* = float
  Power* = float
  Duration* = float
  Interval* = float
  SocTab* = seq[Voltage]

  Tab*[A, B] = seq[(A, B)]


proc Q_from_Ah*(Ah: float): Charge = 
  return Ah * 3600.0  


proc Q_to_Ah*(Q: Charge): float = 
  return Q / 3600.0



proc interpolate*[A, B](tab: Tab[A, B], x: A): B =
  let n = len(tab)
  if x <= tab[0][0]:
    return tab[0][1]
  if x >= tab[n - 1][0]:
    return tab[n - 1][1]
  for i in 0 ..< n - 1:
    if tab[i][0] <= x and x <= tab[i + 1][0]:
      let f1 = tab[i][1]
      let f2 = tab[i + 1][1]
      let x1 = tab[i][0]
      let x2 = tab[i + 1][0]
      return f1 + (f2 - f1) * (x - x1) / (x2 - x1)
  raise newException(ValueError, "Interpolation failed")



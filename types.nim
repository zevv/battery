

type 

  Bfloat* = float32

  Factor* = Bfloat
  Soc* = Bfloat
  Soh* = Bfloat
  Temperature* = Bfloat
  Resistance* = Bfloat
  Capacitance* = Bfloat
  Voltage* = Bfloat
  Current* = Bfloat
  Charge* = Bfloat
  Power* = Bfloat
  Duration* = Bfloat
  Interval* = Bfloat
  SocTab* = seq[Voltage]
  LutFn*[A, B] = proc(x: A): B

  Tab*[A, B] = seq[(A, B)]


proc Q_from_Ah*(Ah: Bfloat): Charge = 
  return Ah * 3600.0  


proc Q_to_Ah*(Q: Charge): Bfloat = 
  return Q / 3600.0




proc interpolate_hermite(tab: Tab[float, float], x: float): float =
  let n = tab.len
  if x <= tab[0][0]: return tab[0][1]
  if x >= tab[n - 1][0]: return tab[n - 1][1]
  var i = 0
  while tab[i+1][0] < x: inc i
  let (x1, y1) = tab[i]
  let (x2, y2) = tab[i+1]
  let dx = x2 - x1
  var m1, m2: float
  if i == 0: m1 = (y2 - y1) / dx
  else: m1 = (y2 - tab[i-1][1]) / (x2 - tab[i-1][0])
  if i >= n - 2: m2 = (y2 - y1) / dx
  else: m2 = (tab[i+2][1] - y1) / (tab[i+2][0] - x1)
  let t = (x - x1) / dx
  let t2 = t * t; let t3 = t2 * t
  let h1 =  2*t3 - 3*t2 + 1
  let h2 = -2*t3 + 3*t2
  let h3 = (t3 - 2*t2 + t) * dx
  let h4 = (t3 -   t2) * dx
  return h1*y1 + h2*y2 + h3*m1 + h4*m2


proc mklut*(tab: Tab[float, float]): LutFn[Bfloat, Bfloat] =
  let x1 = tab[0][0]
  let x2 = tab[^1][0]
  
  const n = 128
  var spline_lut: array[n + 1, Bfloat]

  let step_size = (x2 - x1) / n.Bfloat
  for i in 0..n:
    let current_x = x1 + i.Bfloat * step_size
    spline_lut[i] = interpolate_hermite(tab, current_x)
  
  result = proc(x: Bfloat): Bfloat =
    if x <= x1:
      return spline_lut[0]
    if x >= x2:
      return spline_lut[n]
    let f_idx = (x - x1) / (x2 - x1) * n.Bfloat
    let idx = int(f_idx)
    let frac = f_idx - idx.Bfloat
    return spline_lut[idx] + (spline_lut[idx + 1] - spline_lut[idx]) * frac
    


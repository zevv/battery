import std / [ os, strformat, sequtils, strutils ]

import battery

proc gen_gnuplot*(sim: Simulation, fname: string) =
  let fd = open(fname, fmWrite)
  proc l(s: string) =
    fd.write(s & "\n")

  l("#!/usr/bin/gnuplot -p")
  l("")
  l("reset")
  l("set grid")
  l("set key off")
  l("set multiplot layout 5, 1")
  l("set lmargin at screen 0.08")
  #l("set noxtics")
  l("# 0 margin between multiplots")
  l("set tmargin 1")
  l("set bmargin 1")
  l("set offsets graph 0, 0, 0.05, 0.05")
  l("")

  # Emit inline data blocks for all cells

  for mid, module in sim.pack.modules:
    for cid, cell in module.cells:
      l(&"$cell_{mid:02}_{cid:02} << EOD")
      cell.fd_log.setFilePos(0)
      for line in cell.fd_log.lines:
        l(line)
      l("EOD")
      l("")

  # Emit plotting commands

  proc gen_graph(gs: openArray[tuple[col: int, ylabel: string]], pres: openArray[string]) =
    var ts: seq[string]
    for i, g in gs:
      let lt = if i == 0: "1" else: "2"
      for mid, module in sim.pack.modules:
        for cid, cell in module.cells:
          ts.add(&""" $cell_{mid:02}_{cid:02} u 1:{g.col} w l dt {lt} """ )
    for pre in pres:
      l(pre)
    l(&"""set ylabel "{gs[0].ylabel}"""")
    l(&"""plot {ts.join(", ")}""")
    l("")

  gen_graph([ (2, "I (A)",      )], [ "unset yrange" ])
  gen_graph([ (3, "U (V)",      )], [ "set yrange [2.3:4.4]" ])
  gen_graph([ (4, "SOC (%)",    )], [ "set yrange [-0.1:1.1]" ])
  gen_graph([ (5, "T_core (°C)",),
              (6, "T_case (°C)",)], [ "set yrange [19:30] " ])
  gen_graph([ (7, "SOH (%)",    )], [ "set yrange [-0.1:1.1]" ])

  fd.write("unset multiplot\n")
  fd.close()



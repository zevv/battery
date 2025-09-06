
SRC := $(wildcard *.nim)

battery: $(SRC)
	nim --mm:none -d:release --debugger:native c main.nim

clean:
	rm -f main perf.data perf.data.old battery.gp

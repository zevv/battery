
SRC := $(wildcard *.nim)

battery: $(SRC)
	nim --mm:none -d:release --debugger:native c battery.nim

clean:
	rm -f battery

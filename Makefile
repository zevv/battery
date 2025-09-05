

battery: battery.nim
	nim --mm:none -d:release --debugger:native c battery.nim

clean:
	rm -f battery

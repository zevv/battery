

battery: main.nim battery.nim
	nim --mm:none -d:release --debugger:native c main.nim

clean:
	rm -f battery

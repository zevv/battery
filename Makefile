
SRC := $(wildcard *.nim)

NIMFLAGS += --mm:refc
NIMFLAGS += -d:release
NIMFLAGS += --debugger:native
NIMFLAGS += --passC:-ffast-math
NIMFLAGS += --passC:-march=native

battery: $(SRC) Makefile
	nim $(NIMFLAGS) c main.nim

clean:
	rm -f main perf.data perf.data.old battery.gp

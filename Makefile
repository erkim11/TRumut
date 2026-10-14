CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall
CPPFLAGS += -Isrc
LIBS      = -lz

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
LIBS += -lstdc++fs
endif

PROG = trcall
CMD  = allele prep-catalog build-pon use-pon cancer trio mutual
SRC  = src/main.cpp $(CMD:%=src/cmd/%.cpp) src/lib/catalog.cpp src/lib/input.cpp
HDR  = src/lib/catalog.h src/lib/input.h src/lib/calling.h

all: $(PROG)

$(PROG): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -pthread -o $@ $(SRC) $(LIBS)

clean:
	rm -f $(PROG)

.PHONY: all clean

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
CXXFLAGS += -ffp-contract=off
CPPFLAGS += -I..
LDLIBS   := -lz

SRC := PurpleApplication.cpp PurpleInput.cpp ../common/CpDump.cpp ../common/HumanChromosome.cpp
OBJ := $(SRC:.cpp=.o)

purple_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) purple_port

.PHONY: clean

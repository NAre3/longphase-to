CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDLIBS   := -lz

SRC := AmberApplication.cpp AmberSitesFile.cpp CpDump.cpp HumanChromosome.cpp
OBJ := $(SRC:.cpp=.o)

amber_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) amber_port

.PHONY: clean

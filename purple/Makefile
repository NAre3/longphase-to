CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
CXXFLAGS += -ffp-contract=off -MMD -MP
CPPFLAGS += -I..
LDLIBS   := -lz -pthread

SRC := PurpleApplication.cpp PurpleInput.cpp PurpleSegmentation.cpp PurpleObserved.cpp PurpleFitting.cpp PurpleCopyNumber.cpp PurpleSummary.cpp PurpleWriters.cpp ../common/CpDump.cpp ../common/HumanChromosome.cpp
OBJ := $(SRC:.cpp=.o)
DEP := $(OBJ:.o=.d)

purple_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(DEP) purple_port

-include $(DEP)

.PHONY: clean

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -ffp-contract=off
# 使用 LongPhase-TO in-tree 的 htslib 1.16（行為對照表 §4 引用的即為此版本）
HTSLIB   := ../htslib
CPPFLAGS += -I$(HTSLIB)
LDLIBS   := $(HTSLIB)/libhts.a -lz -lm -lbz2 -llzma -ldeflate -lpthread

SRC := AmberApplication.cpp AmberOutput.cpp AmberSitesFile.cpp BamEvidenceReader.cpp CpDump.cpp \
       HumanChromosome.cpp SamRecordView.cpp TumorFilters.cpp \
       NoiseFloor.cpp CommonsMath.cpp Segmentation.cpp
OBJ := $(SRC:.cpp=.o)

amber_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) amber_port

.PHONY: clean

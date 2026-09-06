CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
# -ffp-contract=off 附加而非併入 CXXFLAGS：環境若帶了自己的 CXXFLAGS，這個旗標不能跟著消失。
# 理由不是與 Java 對齊（已不追求逐位元相同），而是本程式自己的輸出不得隨編譯器與 -march 飄移
# ——PURPLE 直接消費 amber.baf.tsv.gz（4 位小數）與 amber.baf.pcf（#.####）。
# 實測：baseline x86-64 下無作用（ISA 無 FMA），加 -march=native 且允許收縮時 61366 組 CDF 有 5988 組改變。
CXXFLAGS += -ffp-contract=off
# 使用 LongPhase-TO in-tree 的 htslib 1.16（行為對照表 §4 引用的即為此版本）
HTSLIB   := ../htslib
CPPFLAGS += -I$(HTSLIB) -I..
LDLIBS   := $(HTSLIB)/libhts.a -lz -lm -lbz2 -llzma -ldeflate -lpthread

# 共用核心在 ../common（namespace lp），AMBER 與 COBALT 兩邊共用同一份物件檔來源。
COMMON_SRC := ../common/CpDump.cpp ../common/HumanChromosome.cpp ../common/SamRecordView.cpp \
              ../common/Segmentation.cpp

SRC := AmberApplication.cpp AmberOutput.cpp AmberSitesFile.cpp BamEvidenceReader.cpp \
       TumorFilters.cpp NoiseFloor.cpp CommonsMath.cpp Segmentation.cpp
OBJ := $(SRC:.cpp=.o) $(COMMON_SRC:.cpp=.o)

amber_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) amber_port

.PHONY: clean

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
# -ffp-contract=off 附加而非併入 CXXFLAGS：環境若帶了自己的 CXXFLAGS，這個旗標不能跟著消失。
# 理由不是與 Java 對齊（已不追求逐位元相同，spec v2 §1.1 D1），而是本程式自己的輸出
# 不得隨編譯器與 -march 飄移——否則同一份原始碼在不同 build 下給出不同結果，
# 破壞 spec §6 的 rerun_pass。（AMBER 階段實測：baseline x86-64 無 FMA 故無作用，
# 加 -march=native 且允許收縮時 61366 組 CDF 有 5988 組改變。）
CXXFLAGS += -ffp-contract=off
# 使用 LongPhase-TO in-tree 的 htslib（與 AMBER 端同一份）
HTSLIB   := ../htslib
CPPFLAGS += -I$(HTSLIB) -I..
LDLIBS   := $(HTSLIB)/libhts.a -lz -lm -lbz2 -llzma -ldeflate -lpthread

# 共用核心在 ../common（namespace lp），AMBER 與 COBALT 兩邊共用同一份物件檔來源。
COMMON_SRC := ../common/CpDump.cpp ../common/HumanChromosome.cpp \
              ../common/Segmentation.cpp

SRC := CobaltApplication.cpp GcProfile.cpp ReadDepth.cpp Regions.cpp WindowStatuses.cpp

OBJ := $(SRC:.cpp=.o) $(COMMON_SRC:.cpp=.o)

cobalt_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) cobalt_port

.PHONY: clean

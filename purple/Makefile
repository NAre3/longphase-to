CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
CXXFLAGS += -ffp-contract=off -MMD -MP
CPPFLAGS += -I..
LDLIBS   := -lz -pthread

# 物件檔一律放進本 target 自己的 build/，**不寫回原始碼目錄**。
# 理由：../common/*.cpp 由 amber_port、cobalt_port、purple_port 三邊共用，
# 若各自編到 common/*.o 就會互相覆蓋，且 make 看到 .o 比來源新就不重編——
# 於是某支執行檔究竟用哪組旗語建成，取決於「誰後跑」。三支執行檔是保真度
# 驗證的依據，這種靜默的建置不確定性會直接動搖驗證本身。
# 頂層 Makefile 早已為同一原因改用 build/obj（該檔 EXP-I03 的註解）。
# 另一個好處：make -C purple clean 不再刪掉 amber/cobalt 還在用的物件。
OBJDIR := build

SRC := PurpleApplication.cpp PurplePipeline.cpp PurpleInput.cpp PurpleSegmentation.cpp PurpleObserved.cpp PurpleFitting.cpp PurpleCopyNumber.cpp PurpleSummary.cpp PurpleWriters.cpp
# PurpleInputAdapter.cpp **刻意不在此列**：它相依 amber/ 與 cobalt/ 的標頭（進而 htslib），
# purple_port 既不需要也沒有連結 htslib。只有頂層整合建置會編它。
COMMON_SRC := CpDump.cpp HumanChromosome.cpp

OBJ := $(SRC:%.cpp=$(OBJDIR)/%.o) $(COMMON_SRC:%.cpp=$(OBJDIR)/common/%.o)
DEP := $(OBJ:.o=.d)

purple_port: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/common/%.o: ../common/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) purple_port

-include $(DEP)

.PHONY: clean

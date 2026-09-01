// 開發用 harness：直接載入 CP-A5.tsv 跑 noise floor 這一段，不重掃 BAM。
//
// **這不是移植的一部分**，與 -debug_only_chr 同樣是縮短迭代週期的工具：
// CP-A5 的內容已於 EXP-005 逐位元組驗證，且它正是 TumorOnlyPurityAnalysis 的輸入
// （AmberApplication.java:270 的 rawData）。正式記錄的執行仍為完整流程。
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "../NoiseFloor.h"

int main(int argc, char **argv)
{
    if(argc < 2){
        std::fprintf(stderr, "usage: noisefloor_from_cp_a5 <CP-A5.tsv> [grid_scores_out.tsv]\n");
        return 2;
    }

    std::FILE *in = std::fopen(argv[1], "r");
    if(in == nullptr){ std::perror("fopen"); return 2; }

    char line[1024];
    if(std::fgets(line, sizeof(line), in) == nullptr){ return 2; }

    std::vector<amber::PositionEvidence> evidence;
    evidence.reserve(800000);

    while(std::fgets(line, sizeof(line), in) != nullptr){
        char chr[64], ref[8], alt[8];
        int pos, idx, depth, indel, refSup, altSup, bqf, mqf, stf;
        if(std::sscanf(line, "%63[^\t]\t%d\t%d\t%7[^\t]\t%7[^\t]\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                chr, &pos, &idx, ref, alt, &depth, &indel, &refSup, &altSup, &bqf, &mqf, &stf) != 12){
            continue;
        }
        amber::PositionEvidence pe;
        pe.chromosome = chr;
        pe.position = pos;
        pe.ref = ref[0];
        pe.alt = alt[0];
        pe.readDepth = depth;
        pe.indelCount = indel;
        pe.refSupport = refSup;
        pe.altSupport = altSup;
        pe.baseQualFiltered = bqf;
        pe.mapQualFiltered = mqf;
        pe.seqTechFiltered = stf;
        evidence.push_back(std::move(pe));
    }
    std::fclose(in);

    const amber::NoiseFloorResult r = amber::computeNoiseFloor(evidence);

    std::printf("evidencePoints\t%zu\n", r.evidencePoints);
    std::printf("afterImmuneFilter\t%zu\n", r.evidencePointsAfterImmuneFilter);
    std::printf("baselineHetGnomadFrequency\t%.3f\n", r.baselineHetGnomadFrequency);
    for(const amber::PeakDiagnostic &d : r.maximaDiagnostics){
        std::printf("peak\t%.3f\tscore=%.3f\thomProportion=%.3f\tchrArmAUC=%.3f\tmutationAUC=%.3f\tcaptured=%d\t%s\n",
                d.vaf, d.score, d.homozygousProportion, d.chrArmAuc, d.mutationAuc,
                d.capturedPoints, d.classification.c_str());
    }
    std::printf("contamination\t%.17g\n", r.contamination);
    std::printf("noiseFloor\t%.17g\n", r.noiseFloor);

    if(argc >= 3){
        std::FILE *out = std::fopen(argv[2], "w");
        if(out != nullptr){
            std::fprintf(out, "vaf\tscore\n");
            for(const auto &g : r.gridScores){
                std::fprintf(out, "%.3f\t%.17g\n", g.first, g.second);
            }
            std::fclose(out);
        }
    }

    return 0;
}

#ifndef BAF_EVIDENCE_H
#define BAF_EVIDENCE_H

// --baf-source=vcf／vcf-pon：以輸入 VCF 的 het SNP 為 BAF 位點，從 readVariantVec
// （BamParser::get_snp 記下的每條 read 在每個 SNP 的 allele）計數，產出與 AMBER prescan＋
// 共用走訪相同型別的 amber::PositionEvidence，交給不變的 amber::postscan。
//
// 放在頂層而不放在 amber/：amber/ 不得 include longphase-to 的型別（amber_port 要能單獨建置）。
//
// 與 AMBER 的已知差異（刻意不修）：
//   - get_snp 不記錄在該位點有 deletion 或非 REF/ALT 鹼基的 read，readDepth 不含這些 read。
//   - longphase-to 在 get_snp 之前已丟掉 MAPQ < -q（預設 1）的 read，mapQualFiltered 只計
//     MAPQ 介於 -q 與 minMappingQuality-1 的 read。

#include <map>
#include <string>
#include <vector>

#include "Util.h"
#include "ParsingBam.h"
#include "amber/PositionEvidence.h"
#include "common/ChrBaseRegion.h"

enum class BafSiteOrigin { NOT_SOMATIC, PON_ONLY };   // vcf, vcf-pon

// NOT_SOMATIC: origin != SOMATIC.  PON_ONLY: origin == PON.
bool isBafOrigin(VariantOriginType origin, BafSiteOrigin mode);

// "vcf" → NOT_SOMATIC，"vcf-pon" → PON_ONLY。其他值為呼叫端錯誤（Phasing.cpp 已擋）。
BafSiteOrigin bafSiteOriginFromName(const std::string &bafSource);

std::vector<amber::PositionEvidence> evidenceFromReadVariants(
        const std::string &contig,
        const std::map<int, RefAlt> &variants,             // snpFile.getVariants(contig)
        const std::vector<ReadVariant> &reads,             // readVariantVec of this contig
        const std::vector<lp::ChrBaseRegion> &excluded,    // empty if --amber-excluded-bed is not given
        BafSiteOrigin mode,
        int minBaseQuality,                                // --amber-min-base-quality, default 13
        int minMappingQuality);                            // --amber-min-map-quality, default 50

// PURPLE LOH 補回的位點清單（取代 --amber-loci panel）：輸入 VCF 中 SNP、isBafOrigin 為真、
// 不在排除區間內的記錄，het 與 hom 都收（LOH 會把 het 變成 hom）。
// key = stripChrPrefix(contig)，只收 1–22；value = 1-based 位置，已排序。
std::map<std::string, std::vector<int>> lohFillPositionsFromVcf(
        SnpParser &snpFile,
        const std::vector<std::string> &contigs,
        const std::vector<lp::ChrBaseRegion> &excluded,
        BafSiteOrigin mode);

#endif

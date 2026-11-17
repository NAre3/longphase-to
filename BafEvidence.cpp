#include "BafEvidence.h"

#include <algorithm>
#include <cstdlib>
#include <unordered_map>

#include "common/HumanChromosome.h"

namespace {

// 與 amber::prescan 的 CP-A2 同一個判定：逐區間線性掃描，命中即 break（amber/AmberPipeline.cpp）。
bool isExcluded(const std::vector<lp::ChrBaseRegion> &excluded, const std::string &contig, int position1){
    for(const lp::ChrBaseRegion &region : excluded){
        if(region.containsPosition(contig, position1)){
            return true;
        }
    }
    return false;
}

bool isSnp(const RefAlt &refAlt){
    return refAlt.Ref.size() == 1 && refAlt.Alt.size() == 1;
}

bool isAutosome(const std::string &strippedChromosome){
    if(strippedChromosome.empty() || strippedChromosome.size() > 2){
        return false;
    }
    for(char c : strippedChromosome){
        if(c < '0' || c > '9'){
            return false;
        }
    }
    const int n = std::atoi(strippedChromosome.c_str());
    return n >= 1 && n <= 22;
}

}

bool isBafOrigin(VariantOriginType origin, BafSiteOrigin mode){
    if(mode == BafSiteOrigin::NOT_SOMATIC){
        return origin != SOMATIC;
    }
    return origin == PON;
}

BafSiteOrigin bafSiteOriginFromName(const std::string &bafSource){
    return bafSource == "vcf-pon" ? BafSiteOrigin::PON_ONLY : BafSiteOrigin::NOT_SOMATIC;
}

std::vector<amber::PositionEvidence> evidenceFromReadVariants(
        const std::string &contig,
        const std::map<int, RefAlt> &variants,
        const std::vector<ReadVariant> &reads,
        const std::vector<lp::ChrBaseRegion> &excluded,
        BafSiteOrigin mode,
        int minBaseQuality,
        int minMappingQuality)
{
    std::vector<amber::PositionEvidence> sites;
    if(!lp::isHumanChromosome(contig)){
        return sites;
    }

    // ---- 2a：選位點（variants 以 0-based pos 為 key，std::map 已依位置排序）----
    std::unordered_map<int, std::size_t> siteIndex;   // 0-based pos → sites 的索引
    for(const auto &entry : variants){
        const int pos = entry.first;
        const RefAlt &refAlt = entry.second;
        if(!isSnp(refAlt)){
            continue;
        }
        if(refAlt.homozygous){
            continue;
        }
        if(!isBafOrigin(refAlt.originType, mode)){
            continue;
        }
        if(isExcluded(excluded, contig, pos + 1)){
            continue;
        }
        amber::PositionEvidence pe;
        pe.chromosome = contig;
        pe.position = pos + 1;
        pe.ref = refAlt.Ref[0];
        pe.alt = refAlt.Alt[0];
        siteIndex.emplace(pos, sites.size());
        sites.push_back(pe);
    }

    // ---- 2b：計數。每條 read 權重 1（含 fakeRead），不用 processReadVariants 的 0.01 權重 ----
    for(const ReadVariant &read : reads){
        for(const Variant &v : read.variantVec){
            if(v.quality < 0){
                continue;   // INDEL、SV、甲基化等型別；>= 0 為 SNP 的 base quality
            }
            auto found = siteIndex.find(v.position);
            if(found == siteIndex.end()){
                continue;
            }
            amber::PositionEvidence &pe = sites[found->second];
            pe.readDepth += 1;
            const bool mapQualLow = read.mapping_quality < minMappingQuality;
            const bool baseQualLow = v.quality < minBaseQuality;
            // AMBER 兩個 filter 各自計數，一條 read 可兩者都加。
            if(mapQualLow){
                pe.mapQualFiltered += 1;
            }
            if(baseQualLow){
                pe.baseQualFiltered += 1;
            }
            if(mapQualLow || baseQualLow){
                continue;
            }
            if(v.allele == REF_ALLELE){
                pe.refSupport += 1;
            }else if(v.allele == ALT_ALLELE){
                pe.altSupport += 1;
            }
        }
    }

    // indelCount、seqTechFiltered 保持 0；readDepth == 0 的位點保留，由 amber::postscan 移除。
    // sites 依 std::map 的走訪順序建立，已依位置排序。
    return sites;
}

std::map<std::string, std::vector<int>> lohFillPositionsFromVcf(
        SnpParser &snpFile,
        const std::vector<std::string> &contigs,
        const std::vector<lp::ChrBaseRegion> &excluded,
        BafSiteOrigin mode)
{
    std::map<std::string, std::vector<int>> result;
    for(const std::string &contig : contigs){
        const std::string chromosome = lp::stripChrPrefix(contig);
        if(!isAutosome(chromosome)){
            continue;
        }
        const std::map<int, RefAlt> *variants = snpFile.getVariants(contig);
        if(variants == nullptr){
            continue;
        }
        for(const auto &entry : *variants){
            const RefAlt &refAlt = entry.second;
            if(!isSnp(refAlt) || !isBafOrigin(refAlt.originType, mode)){
                continue;
            }
            if(isExcluded(excluded, contig, entry.first + 1)){
                continue;
            }
            result[chromosome].push_back(entry.first + 1);
        }
    }
    for(auto &entry : result){
        std::sort(entry.second.begin(), entry.second.end());
    }
    return result;
}

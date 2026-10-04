#include "Regions.h"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>   // KS_SEP_LINE

#include "CobaltConstants.h"
#include "../common/HumanChromosome.h"
#include "../common/Segmentation.h"

namespace cobalt {

namespace {

std::vector<std::string> splitTab(const std::string &line)
{
    std::vector<std::string> out;
    std::size_t start = 0;
    while(true)
    {
        std::size_t hit = line.find('\t', start);
        if(hit == std::string::npos){ out.push_back(line.substr(start)); break; }
        out.push_back(line.substr(start, hit - start));
        start = hit + 1;
    }
    return out;
}

}

std::vector<ExcludedRegion> loadExcludedRegions(const std::string &path)
{
    htsFile *fp = hts_open(path.c_str(), "r");
    if(fp == nullptr){ throw std::runtime_error("cannot open excluded regions: " + path); }

    std::vector<ExcludedRegion> out;
    kstring_t line = KS_INITIALIZE;
    bool first = true;
    while(hts_getline(fp, KS_SEP_LINE, &line) >= 0)
    {
        if(first){ first = false; continue; }              // ExcludedRegionsFile：跳過表頭
        if(line.l == 0){ continue; }
        while(line.l > 0 && line.s[line.l - 1] == '\r'){ line.s[--line.l] = '\0'; }
        std::vector<std::string> f = splitTab(std::string(line.s, line.l));
        if(f.size() < 3){ continue; }
        ExcludedRegion r;
        r.chromosome = f[0];                                // 不 trim：Java 端也沒有
        r.start = std::atoi(f[1].c_str());
        r.end = std::atoi(f[2].c_str());
        out.push_back(r);
    }
    ks_free(&line);
    hts_close(fp);
    return out;
}

std::size_t DiploidRegions::totalEntries() const
{
    std::size_t n = 0;
    for(const auto &v : byChromosome){ n += v.size(); }
    return n;
}

const std::vector<DiploidStatus> *DiploidRegions::find(const std::string &shortName) const
{
    for(std::size_t i = 0; i < chromosomes.size(); ++i)
    {
        if(chromosomes[i] == shortName){ return &byChromosome[i]; }
    }
    return nullptr;
}

DiploidRegions loadDiploidRegions(const std::string &bedPath)
{
    htsFile *fp = hts_open(bedPath.c_str(), "r");
    if(fp == nullptr){ throw std::runtime_error("cannot open diploid bed: " + bedPath); }

    std::vector<std::string> keys;
    std::vector<std::vector<DiploidStatus>> lists;

    // DiploidRegionLoader 的狀態機（第 55–86 行）
    std::string currentContig;
    bool haveContig = false;
    int position = 0;
    std::size_t currentIdx = 0;

    kstring_t line = KS_INITIALIZE;
    while(hts_getline(fp, KS_SEP_LINE, &line) >= 0)
    {
        if(line.l == 0){ continue; }
        while(line.l > 0 && line.s[line.l - 1] == '\r'){ line.s[--line.l] = '\0'; }
        std::string s(line.s, line.l);
        if(s.rfind("track", 0) == 0 || s.rfind("browser", 0) == 0 || s[0] == '#'){ continue; }
        std::vector<std::string> f = splitTab(s);
        if(f.size() < 3){ continue; }

        const std::string contig = f[0];
        // htsjdk BEDCodec：getStart() = bedStart + 1（1-based 含端點）、getEnd() = bedEnd
        const int bedStart = std::atoi(f[1].c_str()) + 1;
        const int bedEnd = std::atoi(f[2].c_str());

        if(!haveContig || contig != currentContig)
        {
            currentContig = contig;
            haveContig = true;
            position = 1;                                    // 換染色體重設為 1，不是 0
            const std::string shortName = lp::stripChrPrefix(contig);
            currentIdx = keys.size();
            for(std::size_t i = 0; i < keys.size(); ++i)
            {
                if(keys[i] == shortName){ currentIdx = i; break; }
            }
            if(currentIdx == keys.size()){ keys.push_back(shortName); lists.emplace_back(); }
        }

        std::vector<DiploidStatus> &dst = lists[currentIdx];
        // createRatio：先補非二倍體，再放二倍體。兩個迴圈皆為嚴格 <
        for(int p = position; p < bedStart; p += WINDOW_SIZE)
        {
            dst.push_back(DiploidStatus{contig, p, p + WINDOW_SIZE - 1, false});
        }
        for(int q = bedStart; q < bedEnd; q += WINDOW_SIZE)
        {
            dst.push_back(DiploidStatus{contig, q, q + WINDOW_SIZE - 1, true});
        }
        position = bedEnd + 1;
    }
    ks_free(&line);
    hts_close(fp);

    std::vector<std::size_t> order(keys.size());
    for(std::size_t i = 0; i < order.size(); ++i){ order[i] = i; }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b){
        return lp::chromosomeOrdinal(keys[a]) < lp::chromosomeOrdinal(keys[b]);
    });

    DiploidRegions out;
    for(std::size_t i : order)
    {
        out.chromosomes.push_back(keys[i]);
        out.byChromosome.push_back(std::move(lists[i]));
    }
    return out;
}

}

#include "GcProfile.h"

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

// String.trim()：移除 <= U+0020 的前後字元。本檔只會出現空白與 \r，行為等同。
std::string trim(const std::string &s)
{
    std::size_t b = 0, e = s.size();
    while(b < e && static_cast<unsigned char>(s[b]) <= ' '){ ++b; }
    while(e > b && static_cast<unsigned char>(s[e - 1]) <= ' '){ --e; }
    return s.substr(b, e - b);
}

// String.split("\t")：Java 會丟棄**尾端**的空字串，但本檔每行五欄且無尾端空欄，
// 兩者在此等價。仍保留全部欄位以免多欄時行為分歧。
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

bool GcProfile::isMappable() const
{
    return doublesGreaterOrEqual(mappablePercentage, GC_MAPPABLE_THRESHOLD);
}

const std::vector<GcProfile> *GcProfileData::find(const std::string &shortName) const
{
    for(std::size_t i = 0; i < chromosomes.size(); ++i)
    {
        if(chromosomes[i] == shortName){ return &byChromosome[i]; }
    }
    return nullptr;
}

GcProfileData loadGcProfile(const std::string &path)
{
    htsFile *fp = hts_open(path.c_str(), "r");
    if(fp == nullptr){ throw std::runtime_error("cannot open gc profile: " + path); }

    // 以短名為鍵累積，最後依 ordinal 排序（對應 dump 端的 keys.sort(comparingInt(ordinal))）
    std::vector<std::string> keys;
    std::vector<std::vector<GcProfile>> lists;

    kstring_t line = KS_INITIALIZE;
    while(hts_getline(fp, KS_SEP_LINE, &line) >= 0)
    {
        if(line.l == 0){ continue; }
        while(line.l > 0 && line.s[line.l - 1] == '\r'){ line.s[--line.l] = '\0'; }
        std::vector<std::string> f = splitTab(std::string(line.s, line.l));
        if(f.size() < 5){ continue; }

        GcProfile g;
        g.chromosome = trim(f[0]);
        const int pos = std::atoi(trim(f[1]).c_str());
        g.gcContent          = std::strtod(trim(f[2]).c_str(), nullptr);
        g.nonNPercentage     = std::strtod(trim(f[3]).c_str(), nullptr);
        g.mappablePercentage = std::strtod(trim(f[4]).c_str(), nullptr);
        g.start = pos + 1;
        g.end   = pos + WINDOW_SIZE;

        // loadGCContent：只有 HumanChromosome.contains 為真才放進 multimap
        if(!lp::isHumanChromosome(g.chromosome)){ continue; }

        const std::string shortName = lp::stripChrPrefix(g.chromosome);
        std::size_t idx = keys.size();
        for(std::size_t i = 0; i < keys.size(); ++i)
        {
            if(keys[i] == shortName){ idx = i; break; }
        }
        if(idx == keys.size()){ keys.push_back(shortName); lists.emplace_back(); }
        lists[idx].push_back(g);
    }
    ks_free(&line);
    hts_close(fp);

    std::vector<std::size_t> order(keys.size());
    for(std::size_t i = 0; i < order.size(); ++i){ order[i] = i; }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b){
        return lp::chromosomeOrdinal(keys[a]) < lp::chromosomeOrdinal(keys[b]);
    });

    GcProfileData data;
    data.chromosomes.reserve(order.size());
    data.byChromosome.reserve(order.size());
    for(std::size_t i : order)
    {
        data.chromosomes.push_back(keys[i]);
        data.byChromosome.push_back(std::move(lists[i]));
    }
    return data;
}

}

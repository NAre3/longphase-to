#include "PurplePhasedBaf.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include "common/CpDump.h"
#include "common/HumanChromosome.h"

namespace purple {

namespace {

constexpr int CELL_SIZE = 100000;
constexpr double MIN_CELL_DEPTH = 20;
const double KAPPAS[] = {1, 2, 4, 8, 12, 16, 24, 32};
const double SWITCHES[] = {0.001, 0.003, 0.01, 0.03, 0.1};

// FGRID = np.round(np.arange(0.5, 0.9951, 0.0025), 4)：np.arange 的第 k 個元素是 0.5 + k·0.0025，
// np.round 為 rint(x·1e4)/1e4（round-half-even）。
std::vector<double> buildFGrid(){
    std::vector<double> grid;
    for(int k = 0; 0.5 + k * 0.0025 < 0.9951; ++k){ grid.push_back(std::nearbyint((0.5 + k * 0.0025) * 1e4) / 1e4); }
    return grid;
}

// 與 numpy.logaddexp 相同的算法。
double logAddExp(double x, double y){
    if(x == y){ return x + std::log(2.0); }
    const double t = x - y;
    return t > 0 ? x + std::log1p(std::exp(-t)) : y + std::log1p(std::exp(t));
}

struct Cells
{
    std::vector<double> h;     // hap1 比例
    std::vector<double> n;     // read 數（同一條 read 會在多個位點重複計入）
    std::vector<char> samePs;  // samePs[c] 表示第 c 格與第 c − 1 格同一 PS（c ≥ 1）
};

// 段 BAF 對每個 f 的邊際對數概似。
std::vector<double> segmentLogLik(const Cells &cells, const std::vector<double> &grid, double kappa, double s){
    const size_t m = grid.size(), k = cells.h.size();
    const double ls = std::log(s), lk = std::log(1 - s), lh = std::log(0.5);
    std::vector<double> out(m);
    std::vector<double> lp(k), lm(k);
    for(size_t j = 0; j < m; ++j){
        const double f = grid[j];
        for(size_t c = 0; c < k; ++c){
            const double var = f * (1 - f) * kappa / cells.n[c];
            const double norm = -0.5 * std::log(2 * M_PI * var);
            const double dp = cells.h[c] - f, dm = cells.h[c] - (1 - f);
            lp[c] = norm - dp * dp / (2 * var);
            lm[c] = norm - dm * dm / (2 * var);
        }
        double ap = lh + lp[0], am = lh + lm[0];
        for(size_t c = 1; c < k; ++c){
            double np_, nm_;
            if(cells.samePs[c]){
                np_ = logAddExp(ap + lk, am + ls);
                nm_ = logAddExp(am + lk, ap + ls);
            }else{
                np_ = nm_ = logAddExp(ap, am) + lh;
            }
            ap = np_ + lp[c];
            am = nm_ + lm[c];
        }
        out[j] = logAddExp(ap, am);
    }
    return out;
}

}

PhasedBafResult applyPhasedBaf(std::vector<ObservedRegion *> &regions, const std::vector<AmberBaf> &bafs,
        const PhasedSnvMap &phased, const std::set<std::pair<std::string, int>> &lohFilledKeys)
{
    // AMBER 位點依染色體分組、依位置排序。
    std::map<std::string, std::vector<const AmberBaf *>> sites;
    for(const AmberBaf &b : bafs){ sites[lp::stripChrPrefix(b.chromosome)].push_back(&b); }
    for(auto &entry : sites){
        std::stable_sort(entry.second.begin(), entry.second.end(),
                [](const AmberBaf *a, const AmberBaf *b){ return a->position < b->position; });
    }

    // 每段的小格。
    std::vector<Cells> cells(regions.size());
    std::vector<char> hasCells(regions.size(), 0);
    PhasedBafResult result;
    for(size_t i = 0; i < regions.size(); ++i){
        const ObservedRegion &region = *regions[i];
        const std::string chromosome = lp::stripChrPrefix(region.segment.chromosome);
        result.originalBaf.push_back(region.observedBaf);
        result.lohFilled.push_back(lohFilledKeys.count({chromosome, region.segment.start}) ? 1 : 0);
        const auto siteIter = sites.find(chromosome);
        const auto phaseIter = phased.find(chromosome);
        if(siteIter == sites.end() || phaseIter == phased.end()){ continue; }
        // key (PS, 位置 / 100 kb) → (第一個位點位置, Σ h·d, Σ d)
        std::map<std::pair<int, int>, std::tuple<int, double, double>> acc;
        const auto &v = siteIter->second;
        auto it = std::lower_bound(v.begin(), v.end(), region.segment.start,
                [](const AmberBaf *a, int pos){ return a->position < pos; });
        for(; it != v.end() && (*it)->position <= region.segment.end; ++it){
            const AmberBaf &b = **it;
            if(b.tumorDepth <= 0){ continue; }
            const auto q = phaseIter->second.find(b.position);
            if(q == phaseIter->second.end()){ continue; }
            const double h = q->second.first ? b.tumorBaf : 1 - b.tumorBaf;
            const auto key = std::make_pair(q->second.second, b.position / CELL_SIZE);
            auto found = acc.find(key);
            if(found == acc.end()){ found = acc.emplace(key, std::make_tuple(b.position, 0.0, 0.0)).first; }
            std::get<1>(found->second) += h * b.tumorDepth;
            std::get<2>(found->second) += b.tumorDepth;
        }
        // sorted((第一個位點, PS, H, n))：第一個位點各格不同，所以只依位置排序。
        std::vector<std::tuple<int, int, double, double>> list;
        for(const auto &entry : acc){
            const double n = std::get<2>(entry.second);
            if(n >= MIN_CELL_DEPTH){
                list.emplace_back(std::get<0>(entry.second), entry.first.first, std::get<1>(entry.second) / n, n);
            }
        }
        if(list.empty()){ continue; }
        std::sort(list.begin(), list.end());
        Cells &c = cells[i];
        for(size_t j = 0; j < list.size(); ++j){
            c.h.push_back(std::get<2>(list[j]));
            c.n.push_back(std::get<3>(list[j]));
            c.samePs.push_back(j > 0 && std::get<1>(list[j]) == std::get<1>(list[j - 1]) ? 1 : 0);
        }
        hasCells[i] = 1;
    }

    // (κ, s) 格點搜尋：所有段 max_f 對數概似的總和最大者。
    const std::vector<double> grid = buildFGrid();
    bool first = true;
    double bestTotal = 0;
    for(double kappa : KAPPAS){
        for(double s : SWITCHES){
            double total = 0;
            for(size_t i = 0; i < regions.size(); ++i){
                if(!hasCells[i]){ continue; }
                const auto ll = segmentLogLik(cells[i], grid, kappa, s);
                total += *std::max_element(ll.begin(), ll.end());
            }
            if(first || total > bestTotal){
                first = false; bestTotal = total; result.kappa = kappa; result.switchProbability = s;
            }
        }
    }

    for(size_t i = 0; i < regions.size(); ++i){
        result.cellCount.push_back(hasCells[i] ? static_cast<int>(cells[i].h.size()) : 0);
        if(hasCells[i]){ ++result.regionsWithCells; }
        if(hasCells[i] && !result.lohFilled[i]){
            const auto ll = segmentLogLik(cells[i], grid, result.kappa, result.switchProbability);
            regions[i]->observedBaf = grid[std::max_element(ll.begin(), ll.end()) - ll.begin()];
            ++result.regionsUpdated;
        }else{
            result.cellCount.back() = 0;
        }
        result.phasedBaf.push_back(regions[i]->observedBaf);
    }
    return result;
}

void dumpPhasedBaf(const std::vector<ObservedRegion *> &regions, const PhasedBafResult &result)
{
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(regions.size());
    for(size_t i = 0; i < regions.size(); ++i){
        const auto &s = regions[i]->segment;
        rows.push_back({s.chromosome, s.start, s.chromosome + "\t" + std::to_string(s.start) + "\t" + std::to_string(s.end) + "\t" +
                std::to_string(regions[i]->bafCount) + "\t" + lp::CpDump::num(result.originalBaf[i]) + "\t" +
                lp::CpDump::num(result.phasedBaf[i]) + "\t" + std::to_string(result.cellCount[i]) + "\t" +
                (result.lohFilled[i] ? "true" : "false") + "\t" + lp::CpDump::num(result.kappa) + "\t" +
                lp::CpDump::num(result.switchProbability)});
    }
    lp::CpDump::write("CP-P4b-phased-baf", "chromosome\tstart\tend\tbafCount\tobservedBafOriginal\tobservedBafPhased\tcells\tlohFilled\tkappa\tswitch", rows);
}

}

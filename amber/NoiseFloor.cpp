#include "NoiseFloor.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#include "CommonsMath.h"
#include "HumanChromosome.h"

namespace amber {

namespace {

#include "Centromeres38.inc"

// ---- ChrArmLocator（ChrArmLocator.java:11-21）----
// position < centromere → P 臂，否則 Q 臂。以 (chromosome, arm) 為分類鍵。
int centromereOf(const std::string &chromosome)
{
    const std::string stripped = stripChrPrefix(chromosome);
    for(const auto &entry : CENTROMERES_38){
        if(stripped == entry.first){
            return entry.second;
        }
    }
    throw std::runtime_error("no centromere for chromosome " + chromosome);
}

std::string chrArmOf(const PositionEvidence &pe)
{
    const std::string stripped = stripChrPrefix(pe.chromosome);
    return stripped + (pe.position < centromereOf(pe.chromosome) ? "_P" : "_Q");
}

// ---- CanonicalSnvType（CanonicalSnvType.java:20-58）----
// ref 為 A 或 G 時取互補，再依 C_*/T_* 分類。
char complementBase(char base)
{
    switch(base){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

std::string canonicalSnvType(char ref, char alt)
{
    char canonicalRef = ref;
    char canonicalAlt = alt;

    if(ref == 'A' || ref == 'G'){
        canonicalRef = complementBase(ref);
        canonicalAlt = complementBase(alt);
    }

    // Java 端在 ref/alt 相同或為 N 時丟 IllegalArgumentException；
    // 落到 null 的組合（例如互補後仍非 C/T）會讓 HashMap 以 null 為鍵。
    // 這裡以固定字串代表，等價於一個獨立分類。
    std::string result;
    result += canonicalRef;
    result += '_';
    result += canonicalAlt;
    return result;
}

// ---- PositionEvidence 的兩個 vaf（PositionEvidence.java:70-92）----
double vafOf(const PositionEvidence &pe)
{
    if(pe.readDepth == 0){
        return std::nan("");
    }
    return static_cast<double>(pe.altSupport) / (pe.altSupport + pe.refSupport);
}

double symmetricVafOf(const PositionEvidence &pe)
{
    const int vafDepth = pe.altSupport + pe.refSupport;
    if(vafDepth == 0){
        return std::nan("");
    }
    int k = pe.altSupport;
    if(k > vafDepth / 2){
        k = pe.refSupport;
    }
    return static_cast<double>(k) / vafDepth;
}

// ---- SearchGrid（SearchGrid.java:26-39）----
std::vector<std::pair<double, double>> searchValuesAndSteps()
{
    std::vector<std::pair<double, double>> result;
    double step = PEAK_SEARCH_INITIAL_STEP;
    double current = PEAK_SEARCH_START;

    while(current <= PEAK_SEARCH_END + PEAK_SEARCH_OVERSHOOT){
        // Java 的 Math.round(double) 回傳 long，且為 floor(x + 0.5)
        const double currentValue = std::floor(current * 1000.0 + 0.5) / 1000.0;
        current += step;
        step *= PEAK_SEARCH_STEP_RATIO;
        result.emplace_back(currentValue, step);
    }

    return result;
}

// ---- RangeStepPeak.isCaptured（CandidatePeak.java:84-94）----
bool isCaptured(int depth, int count, double vaf, double step, double symmetricVaf)
{
    const double cdf = cm3::binomialCdf(depth, vaf, count);
    const bool capturedByCdf = cdf > LOWER_CDF_BOUND_FOR_CAPTURE && cdf < UPPER_CDF_BOUND_FOR_CAPTURE;
    const bool capturedByStep = std::fabs(symmetricVaf - vaf) < step;
    return capturedByCdf || capturedByStep;
}

// ---- CandidatePeak（CandidatePeak.java）----
class CandidatePeak
{
public:
    CandidatePeak(double level, double stepToNextLevel)
        : mLevel(level), mStep(stepToNextLevel) {}

    // CandidatePeak.java:35-40
    bool hasSufficientDepthForEventDetection(const PositionEvidence &pe) const
    {
        const int n = pe.refSupport + pe.altSupport;
        return cm3::binomialCdf(n, mLevel / 2, 2) < LOWER_CDF_BOUND_FOR_CAPTURE;
    }

    // CandidatePeak.java:42-82
    void test(const PositionEvidence &pe)
    {
        int n = pe.refSupport + pe.altSupport;
        int k = pe.altSupport;
        if(k > n / 2){
            k = n - k;
        }

        const double symmetricVaf = symmetricVafOf(pe);
        const bool hom = isCaptured(n, k, mLevel, mStep, symmetricVaf);
        const bool het = isCaptured(n, k, mLevel / 2, mStep / 2, symmetricVaf);

        if(hom){
            if(het){
                if(symmetricVaf > mLevel * 0.75){
                    mHomozygousBand.push_back(&pe);
                }else{
                    mHeterozygousBand.push_back(&pe);
                }
            }else{
                mHomozygousBand.push_back(&pe);
            }
        }else if(het){
            mHeterozygousBand.push_back(&pe);
        }
    }

    // allCapturedPoints 在 Java 端是兩個 band 的 HashSet 聯集。
    // 每個 evidence 至多被 test 一次且只會進一個 band，故聯集大小即兩者相加。
    std::size_t numberOfCapturedEvidencePoints() const
    {
        return mHomozygousBand.size() + mHeterozygousBand.size();
    }

    std::vector<const PositionEvidence *> allCapturedPoints() const
    {
        std::vector<const PositionEvidence *> all;
        all.reserve(numberOfCapturedEvidencePoints());
        all.insert(all.end(), mHomozygousBand.begin(), mHomozygousBand.end());
        all.insert(all.end(), mHeterozygousBand.begin(), mHeterozygousBand.end());
        return all;
    }

    double homozygousProportion() const
    {
        return static_cast<double>(mHomozygousBand.size()) / numberOfCapturedEvidencePoints();
    }

    double vaf() const { return mLevel; }

private:
    double mLevel;
    double mStep;
    std::vector<const PositionEvidence *> mHomozygousBand;
    std::vector<const PositionEvidence *> mHeterozygousBand;
};

// ---- AucCalculator + CategoryEvidence + CategoryEvidenceIntegral ----
// AucCalculator.java:21-42、CategoryEvidence.java:44-62、CategoryEvidenceIntegral.java:12-39
class AucCalculator
{
public:
    explicit AucCalculator(std::map<std::string, int> baseline)
        : mBaseline(std::move(baseline)) {}

    double calculateAuc(const std::map<std::string, int> &counts) const
    {
        struct Entry
        {
            std::string category;
            int totalPoints;   // baseline
            int evidencePoints; // sample
            double ratio() const
            {
                if(totalPoints == 0){
                    return std::numeric_limits<double>::max(); // Double.MAX_VALUE
                }
                return static_cast<double>(evidencePoints) / totalPoints;
            }
        };

        std::vector<std::string> categories;
        for(const auto &kv : counts){ categories.push_back(kv.first); }
        for(const auto &kv : mBaseline){
            if(counts.find(kv.first) == counts.end()){ categories.push_back(kv.first); }
        }

        std::vector<Entry> entries;
        entries.reserve(categories.size());
        for(const std::string &category : categories){
            const auto b = mBaseline.find(category);
            const auto s = counts.find(category);
            entries.push_back({category,
                    b == mBaseline.end() ? 0 : b->second,
                    s == counts.end() ? 0 : s->second});
        }

        // CategoryEvidence.compareTo：先比 ratio，相同再比 category。
        // HashSet 的走訪順序在此不影響結果，因為隨後即排序（全序）。
        std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b){
            if(a.ratio() != b.ratio()){
                return a.ratio() < b.ratio();
            }
            return a.category < b.category;
        });

        long totalHits = 0;
        long totalPoints = 0;
        for(const Entry &e : entries){
            totalHits += e.evidencePoints;
            totalPoints += e.totalPoints;
        }

        double heightSoFar = 0.0;
        double area = 0.0;
        for(const Entry &e : entries){
            area += e.totalPoints * heightSoFar;
            area += e.evidencePoints * e.totalPoints / 2.0;
            heightSoFar += e.evidencePoints;
        }

        const double possibleMax = 0.5 * totalHits * totalPoints;
        return area / possibleMax;
    }

private:
    std::map<std::string, int> mBaseline;
};

template<typename Classifier>
std::map<std::string, int> getCounts(const std::vector<const PositionEvidence *> &samples, Classifier classifier)
{
    std::map<std::string, int> counts;
    for(const PositionEvidence *pe : samples){
        ++counts[classifier(*pe)];
    }
    return counts;
}

// ---- RegionsFilter（RegionsFilter.java）：免疫區排除 ----
// ImmuneRegions.java 的 V38 清單（HLA 前後各 5000 bp flank）
struct ImmuneRegion { const char *chromosome; int start; int end; };

const ImmuneRegion IMMUNE_REGIONS_V38[] = {
    // HLA（TRANSCRIPT_FLANK = 5000）
    {"chr6", 29942532 - 5000, 29945870 + 5000},
    {"chr6", 31353872 - 5000, 31367067 + 5000},
    {"chr6", 31268749 - 5000, 31272092 + 5000},
    // IG
    {"chr2", 88857161, 90315836},
    {"chr14", 105586437, 106879844},
    {"chr22", 22026076, 22922913},
    // TR
    {"chr7", 38239580, 38367882},
    {"chr7", 142299177, 142812869},
    {"chr14", 21622293, 22552156},
};

// RegionsFilter.java:18-52 的狀態式掃描：依染色體分組、組內排序，
// 走訪輸入時只檢查「第一個 end >= position 的區間」。輸入需依 (染色體, position) 遞增，
// CP-A5 的 rawData 正是如此（EXP-005 已驗證其排序）。
// 刻意不改成「檢查所有區間」——雖然在不重疊的區間下等價，但那是額外假設。
std::vector<const PositionEvidence *> filterOutImmuneRegions(const std::vector<PositionEvidence> &evidence)
{
    std::map<std::string, std::vector<std::pair<int, int>>> byChromosome;
    for(const ImmuneRegion &r : IMMUNE_REGIONS_V38){
        byChromosome[r.chromosome].emplace_back(r.start, r.end);
    }
    for(auto &entry : byChromosome){
        std::sort(entry.second.begin(), entry.second.end());
    }

    std::vector<const PositionEvidence *> result;
    result.reserve(evidence.size());

    std::string currentChromosome;
    bool haveChromosome = false;
    const std::vector<std::pair<int, int>> *currentRegions = nullptr;
    std::size_t regionIndex = 0;

    for(const PositionEvidence &pe : evidence){
        if(!haveChromosome || pe.chromosome != currentChromosome){
            currentChromosome = pe.chromosome;
            haveChromosome = true;
            const auto found = byChromosome.find(pe.chromosome);
            currentRegions = found == byChromosome.end() ? nullptr : &found->second;
            regionIndex = 0;
        }

        if(currentRegions == nullptr){
            result.push_back(&pe);
            continue;
        }

        while(regionIndex < currentRegions->size() && (*currentRegions)[regionIndex].second < pe.position){
            ++regionIndex;
        }

        if(regionIndex >= currentRegions->size()
                || !(pe.position >= (*currentRegions)[regionIndex].first
                        && pe.position <= (*currentRegions)[regionIndex].second)){
            result.push_back(&pe);
        }
    }

    return result;
}

// ---- LocalMaximaFinder（LocalMaximaFinder.java:22-63）----
struct GridResult
{
    CandidatePeak peak;
    double score;
};

std::vector<std::size_t> findLocalMaxima(const std::vector<double> &scores)
{
    std::vector<std::size_t> maxima;

    std::size_t firstNonZero = scores.size();
    for(std::size_t i = 0; i < scores.size(); ++i){
        if(scores[i] != 0.0){ firstNonZero = i; break; }
    }
    if(firstNonZero == scores.size()){
        return maxima;
    }

    enum class Direction { NONE, UP, DOWN, FLAT };
    Direction previousDirection = Direction::NONE;
    Direction currentDirection = Direction::NONE;
    bool hasPrevious = false;
    std::size_t previousIndex = 0;

    for(std::size_t i = 0; i < scores.size(); ++i){
        if(hasPrevious){
            if(scores[i] > scores[previousIndex]){
                currentDirection = Direction::UP;
            }else if(scores[i] < scores[previousIndex]){
                currentDirection = Direction::DOWN;
            }else{
                currentDirection = Direction::FLAT;
            }

            if(previousDirection == Direction::UP || previousDirection == Direction::FLAT){
                if(currentDirection == Direction::DOWN){
                    maxima.push_back(previousIndex);
                }
            }
        }

        previousDirection = currentDirection;
        previousIndex = i;
        hasPrevious = true;
    }

    // maxima.remove(firstNonZero)：List.remove(Object) 移除**第一個相等的元素**。
    // CandidatePeakEvaluationResult 是 record，equals 比 (candidatePeak, score, message)，
    // 其中 CandidatePeak 未覆寫 equals → 物件識別。故等價於「移除第一個就是該筆結果的元素」。
    for(auto it = maxima.begin(); it != maxima.end(); ++it){
        if(*it == firstNonZero){
            maxima.erase(it);
            break;
        }
    }

    return maxima;
}

}

NoiseFloorResult computeNoiseFloor(const std::vector<PositionEvidence> &evidence)
{
    NoiseFloorResult result;
    result.evidencePoints = evidence.size();

    // TumorOnlyPurityAnalysis.java:67 filterOutExcludedRegions
    const std::vector<const PositionEvidence *> filtered = filterOutImmuneRegions(evidence);
    result.evidencePointsAfterImmuneFilter = filtered.size();

    // TumorOnlyPurityAnalysis.java:150-162 getBaselineHetVariants
    std::vector<const PositionEvidence *> hetVariants;
    for(const PositionEvidence *pe : filtered){
        const double vaf = symmetricVafOf(*pe);
        if(vaf > HET_VAF_LOWER_BOUND && vaf < HET_VAF_UPPER_BOUND){
            hetVariants.push_back(pe);
        }
    }

    // gnomad 頻率：AmberGermlineSites.38.tsv.gz 無 Frequency 欄，全部為 0
    // → baseline 平均為 0，且每個 peak 的兩個 band 平均亦為 0，
    //   PeakGnomadFrequenciesChecker 因此退化為「兩個 band 皆非空」（見下方）。
    //
    // **已知限制（2026-09-02 程式碼審查 F-R1）**：此處是把「頻率全為 0」下的
    // 退化行為寫死，而非移植完整檢查。Java 的 checkGnomadFrequencies 有三個子句
    // （jar bytecode 156-217）：
    //     getN()==0 任一為真         → reject
    //     |low.mean  - baseline| > 0.15 → reject
    //     |high.mean - baseline| > 0.15 → reject
    // 本檔只實作第一句，門檻 0.15 不存在於此移植中。前提「頻率全為 0」沒有被
    // 任何程式碼強制（AmberSitesFile 仍會讀 Frequency 欄，只是下游丟棄）。
    //
    // 已知 owner 判定（2026-09-02）：不加防呆拋錯，僅註記。見 AmberSitesFile.cpp
    // 同一編號的註解，與 research/studies/purple-port-amber-fidelity-v1/CODE-REVIEW-fcb15ca.md
    result.baselineHetGnomadFrequency = 0.0;

    const auto chrArmClassifier = [](const PositionEvidence &pe){ return chrArmOf(pe); };
    const auto snvTypeClassifier = [](const PositionEvidence &pe){
        // TumorOnlyPurityAnalysis.java:39-52
        if(vafOf(pe) <= HET_VAF_UPPER_BOUND){
            return canonicalSnvType(pe.ref, pe.alt);
        }
        return canonicalSnvType(pe.alt, pe.ref);
    };

    const AucCalculator chrArmAuc(getCounts(hetVariants, chrArmClassifier));
    const AucCalculator snvTypeAuc(getCounts(hetVariants, snvTypeClassifier));

    // PeakSearch.java:20-47
    const std::vector<std::pair<double, double>> grid = searchValuesAndSteps();

    std::vector<CandidatePeak> peaks;
    std::vector<double> scores;
    peaks.reserve(grid.size());
    scores.reserve(grid.size());

    for(const auto &levelAndStep : grid){
        CandidatePeak peak(levelAndStep.first, levelAndStep.second);

        // CandidatePeakEvaluation.java:23-42
        std::vector<const PositionEvidence *> testable;
        testable.reserve(filtered.size());
        for(const PositionEvidence *pe : filtered){
            if(peak.hasSufficientDepthForEventDetection(*pe)){
                testable.push_back(pe);
            }
        }

        double score;
        if(testable.size() < static_cast<std::size_t>(MINIMUM_CAPTURED_POINTS)){
            score = 0.0;
        }else{
            for(const PositionEvidence *pe : testable){
                peak.test(*pe);
            }
            if(peak.numberOfCapturedEvidencePoints() < static_cast<std::size_t>(MINIMUM_CAPTURED_POINTS)){
                score = 0.0;
            }else{
                score = static_cast<double>(peak.numberOfCapturedEvidencePoints()) / testable.size();
            }
        }

        peaks.push_back(std::move(peak));
        scores.push_back(score);
        result.gridScores.emplace_back(levelAndStep.first, score);
    }

    const std::vector<std::size_t> maxima = findLocalMaxima(scores);

    // TumorOnlyPurityAnalysis.java:82-111
    std::vector<double> copyNumberPeakVafs;

    for(std::size_t index : maxima){
        const CandidatePeak &peak = peaks[index];
        const std::vector<const PositionEvidence *> captured = peak.allCapturedPoints();

        PeakDiagnostic diagnostic;
        diagnostic.vaf = peak.vaf();
        diagnostic.score = scores[index];
        diagnostic.capturedPoints = static_cast<int>(captured.size());

        // PeakGnomadFrequenciesChecker.java:21-51
        // 所有頻率為 0 → 兩個平均皆為 0、expectedMean 為 0 → 差值 0，必然通過門檻；
        // 唯一可能失敗的是某一側為空（getN() == 0）。
        std::size_t lowCount = 0;
        std::size_t highCount = 0;
        for(const PositionEvidence *pe : captured){
            if(vafOf(*pe) < 0.5){ ++lowCount; }else{ ++highCount; }
        }
        if(lowCount == 0 || highCount == 0){
            diagnostic.homozygousProportion = peak.homozygousProportion();
            diagnostic.chrArmAuc = 0.0;
            diagnostic.mutationAuc = 0.0;
            diagnostic.classification = "gnomad-rejected";
            result.maximaDiagnostics.push_back(diagnostic);
            continue;
        }

        const double homozygousProportion = peak.homozygousProportion();
        const double armAuc = chrArmAuc.calculateAuc(getCounts(captured, chrArmClassifier));
        const double mutAuc = snvTypeAuc.calculateAuc(getCounts(captured, snvTypeClassifier));

        diagnostic.homozygousProportion = homozygousProportion;
        diagnostic.chrArmAuc = armAuc;
        diagnostic.mutationAuc = mutAuc;
        diagnostic.classification = "none";

        // Range.open(lower, upper).contains：開區間，兩端不含
        const bool inHomPropRange = homozygousProportion > HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION
                && homozygousProportion < HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION;

        if(inHomPropRange){
            if(mutAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION
                    && armAuc > CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION){
                diagnostic.classification = "contamination";
                result.contaminationPeakVafs.push_back(peak.vaf());
                result.maximaDiagnostics.push_back(diagnostic);
                continue;
            }
        }

        if(armAuc < CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS
                || (mutAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION
                        && armAuc < CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION)){
            diagnostic.classification = "copy-number";
            copyNumberPeakVafs.push_back(peak.vaf());
        }

        result.maximaDiagnostics.push_back(diagnostic);
    }

    // TumorOnlyPurityAnalysis.cutoff()（同檔 130-138）
    double minCopyNumberPeak = 0.0;
    if(!copyNumberPeakVafs.empty()){
        minCopyNumberPeak = *std::min_element(copyNumberPeakVafs.begin(), copyNumberPeakVafs.end());
    }

    if(doublesIsZero(minCopyNumberPeak)){
        result.noiseFloor = MIN_CUTOFF;
    }else{
        result.noiseFloor = std::fmin(MIN_CUTOFF, minCopyNumberPeak / 3.0);
    }

    // AmberApplication.java:281：contamination 取汙染 peak 的最大 vaf，無則 0.0
    result.contamination = 0.0;
    for(double vaf : result.contaminationPeakVafs){
        result.contamination = std::fmax(result.contamination, vaf);
    }

    return result;
}

}

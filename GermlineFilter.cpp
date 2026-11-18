#include "GermlineFilter.h"

#include <algorithm>
#include <cerrno>
#include <memory>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "purple/PurpleWriters.h"

const char *const kGermlineFilterTag = "LP_germlineHG";

namespace germline_filter {
namespace {

// ---------------------------------------------------------------- 常數（E16 README §0.4 逐項列出來源）
constexpr int kExcludeFlags = 0xF04;      // unmapped／secondary／QC fail／duplicate（E03 §0.4）
constexpr int kSupplementaryFlag = 0x800; // supplementary（E03 §0.4）
constexpr int kMinMapq = 5;               // E03 §0.4
constexpr int kMinBaseQuality = 10;       // E03 §0.4（anchor 與目標位點）
constexpr double kAlpha = 0.05;           // 可判（E03 §0.5、E06 §0.5）
constexpr double kAlphaHom = 1e-3;        // 對照剔除、E07 主設計（E06 §0.5、E09 §0.2）
constexpr double kEpsExclude = 0.001;     // 偽 het 判定的 ε（E08 §0.3、E09 §0.2）
constexpr double kErrFloor = 0.001;       // err、ε 的下限（E03 §0.5、E06 §0.11）
constexpr int kMinCtrlPerSeg = 100;       // 段型對照少於此數用全組合併值（E03 §0.5）
constexpr int kMinGnomadAc = 2;           // gnomAD AC>=2（E06 §0.6）

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kNegInf = -std::numeric_limits<double>::infinity();

// ---------------------------------------------------------------- 與 Python／numpy／scipy 一致的數值函式
// Python round(x)（half-even）與 np.round(x)：預設捨入模式下的 nearbyint
double roundHalfEven(double x){ return std::nearbyint(x); }
// Python round(x, 6)：十進位正確捨入（glibc printf 同樣是精確值上的 half-even）
double pyRound6(double x){
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.6f", x);
    return std::strtod(buffer, nullptr);
}
// np.round(x, 6)：rint(x * 1e6) / 1e6
double npRound6(double x){ return std::nearbyint(x * 1e6) / 1e6; }

double xlogy(double x, double y){ return x == 0 ? 0.0 : x * std::log(y); }
double xlog1py(double x, double y){ return x == 0 ? 0.0 : x * std::log1p(y); }
// scipy.stats.binom._logpmf：gammaln 組合數 + xlogy(k, p) + xlog1py(n-k, -p)
double binomLogpmf(int k, int n, double p){
    if(k < 0 || k > n){ return kNegInf; }
    const double combiln = std::lgamma(n + 1.0) - (std::lgamma(k + 1.0) + std::lgamma(n - k + 1.0));
    return combiln + xlogy(k, p) + xlog1py(n - k, -p);
}
double binomPmf(int k, int n, double p){ return std::exp(binomLogpmf(k, n, p)); }
// P(Bin(n, p) >= r)，即 scipy binom.sf(r - 1, n, p)
double binomTailGe(int r, int n, double p){
    if(r <= 0){ return 1.0; }
    if(r > n){ return 0.0; }
    double s = 0;
    for(int j = r; j <= n; ++j){ s += binomPmf(j, n, p); }
    return std::min(s, 1.0);
}
// scipy logsumexp([a, b], b=.5)
double logHalfSumExp(double a, double b){
    const double m = std::max(a, b);
    if(m == kNegInf){ return kNegInf; }
    return std::log(0.5 * std::exp(a - m) + 0.5 * std::exp(b - m)) + m;
}
double clip(double x, double lo, double hi){ return std::min(std::max(x, lo), hi); }

// ---------------------------------------------------------------- 文字欄位
void splitTabs(const std::string &line, std::vector<std::string> &out){
    out.clear();
    size_t start = 0;
    size_t end = line.size();
    if(end > 0 && line[end - 1] == '\n'){ --end; }
    while(true){
        size_t tab = line.find('\t', start);
        if(tab == std::string::npos || tab >= end){
            out.emplace_back(line, start, end - start);
            break;
        }
        out.emplace_back(line, start, tab - start);
        start = tab + 1;
    }
}
std::vector<std::string> splitChar(const std::string &s, char c){
    std::vector<std::string> out;
    size_t start = 0;
    while(true){
        size_t p = s.find(c, start);
        if(p == std::string::npos){ out.push_back(s.substr(start)); break; }
        out.push_back(s.substr(start, p - start));
        start = p + 1;
    }
    return out;
}
std::string upper(std::string s){
    for(char &c : s){ c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
    return s;
}
// pandas.to_numeric(errors='coerce') 的整數子集：整個字串是整數才算
bool parseInt(const std::string &s, long &out){
    if(s.empty()){ return false; }
    char *end = nullptr;
    errno = 0;
    long v = std::strtol(s.c_str(), &end, 10);
    if(errno != 0 || end == s.c_str() || *end != '\0'){ return false; }
    out = v;
    return true;
}
int baseIndex(char c){
    switch(c){ case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
bool isAutosome(const std::string &chrom){
    std::string c = chrom.rfind("chr", 0) == 0 ? chrom.substr(3) : chrom;
    long v;
    return parseInt(c, v) && v >= 1 && v <= 22 && std::to_string(v) == c;
}

// ---------------------------------------------------------------- PURPLE 段
struct Segment { int start, end; double major, minor; };
struct SegmentIndex {
    std::map<std::string, std::vector<Segment>> byChrom;
    // E03 annotate_cnv：start <= pos <= end 的最後一段（start 排序後 searchsorted right − 1）
    bool find(const std::string &chrom, int pos1, double &mj, double &mn) const {
        mj = mn = kNaN;
        auto it = byChrom.find(chrom);
        if(it == byChrom.end() || it->second.empty()){ return false; }
        const auto &v = it->second;
        auto up = std::upper_bound(v.begin(), v.end(), pos1, [](int p, const Segment &s){ return p < s.start; });
        if(up == v.begin()){ return false; }
        --up;
        if(pos1 > up->end){ return false; }
        mj = up->major; mn = up->minor;
        return true;
    }
};
enum SegType { SEG_LOH = 0, SEG_BAL = 1, SEG_GAIN = 2, SEG_NA = 3 };
const char *segName(int s){ static const char *n[] = {"LOH", "bal", "gain", "NA"}; return n[s]; }
// E03 segtype
int segType(double mj, double mn){
    if(std::isnan(mn)){ return SEG_NA; }
    if(mn < .5){ return SEG_LOH; }
    return roundHalfEven(mj) > roundHalfEven(mn) ? SEG_GAIN : SEG_BAL;
}

// ---------------------------------------------------------------- 位點
struct Anchor {
    int pos0;
    int ref, alt;   // 0..3，非 ACGT 為 -1
    int hapAlt;     // GT 0|1 → 2；1|0 → 1
    long ps;        // 輸出 VCF 的 PS 字串值
    bool adKnown = false;
    long adRef = 0, adN = 0;
    bool baseKnown = false;
    long baseTotal = 0, baseThird = 0;
    bool excluded = false;
    bool atPass = false;
};
struct Target {
    int pos0;
    std::string ref, alt;   // 大寫
    int refB, altB;
    bool isCtrl;
    // 計數：h1r h1a h2r h2a unr una cfr cfa oth vote_sum n_voted
    long cnt[11] = {0};
    long tps = -1;
    int dist = -1;
};
struct ChromSites {
    std::vector<Anchor> anchors;     // 輸出順序
    std::set<int> passPos0;          // FILTER=PASS 的位置（不分類型）
    std::vector<Target> targets;     // PASS SNV ＋ 對照
    std::vector<Target> gtSelected;  // GT 挑到的 SNV（TDF > 0.9）
};

// ---------------------------------------------------------------- pre-pass：從模擬的輸出行取位點
void collectSites(SnpParser &snpFile, const std::string &vcf, ChrPhasingResult &chrPhasingResult,
                  bool highPurity, std::map<std::string, ChromSites> &sites, long &nRecords){
    std::vector<std::string> outF, inF;
    snpFile.scanOutputRecords(vcf, chrPhasingResult, [&](const std::string &input, const std::string &outLine){
        ++nRecords;
        splitTabs(outLine, outF);
        if(outF.size() < 10){ return; }
        const std::string &chrom = outF[0];
        const int pos1 = std::atoi(outF[1].c_str());
        const std::string &flt = outF[6];
        const bool snv = outF[3].size() == 1 && outF[4].size() == 1;
        ChromSites &cs = sites[chrom];
        if(flt == "PASS"){ cs.passPos0.insert(pos1 - 1); }
        if(!snv){ return; }
        if(highPurity){
            if(lineMatchesSomaticGt(outLine)){
                Target t;
                t.pos0 = pos1 - 1; t.ref = upper(outF[3]); t.alt = upper(outF[4]);
                t.refB = baseIndex(t.ref[0]); t.altB = baseIndex(t.alt[0]); t.isCtrl = false;
                cs.gtSelected.push_back(t);
            }
            return;
        }
        if(flt == "PASS"){
            Target t;
            t.pos0 = pos1 - 1; t.ref = upper(outF[3]); t.alt = upper(outF[4]);
            t.refB = baseIndex(t.ref[0]); t.altB = baseIndex(t.alt[0]); t.isCtrl = false;
            cs.targets.push_back(t);
            return;
        }
        if(flt != "LP_PON"){ return; }
        // E03 read_anchor_vcf：GT＝第一個 FORMAT 值，PS 依 FORMAT 位置
        std::vector<std::string> fmt = splitChar(outF[8], ':');
        std::vector<std::string> val = splitChar(outF[9], ':');
        const std::string &gt = val[0];
        if(gt != "0|1" && gt != "1|0"){ return; }
        auto psIt = std::find(fmt.begin(), fmt.end(), std::string("PS"));
        if(psIt == fmt.end()){ return; }
        const size_t psIdx = psIt - fmt.begin();
        if(psIdx >= val.size() || val[psIdx] == "."){ return; }
        long ps;
        if(!parseInt(val[psIdx], ps)){ return; }
        Anchor a;
        a.pos0 = pos1 - 1;
        a.ref = baseIndex(upper(outF[3])[0]);
        a.alt = baseIndex(upper(outF[4])[0]);
        a.hapAlt = gt == "0|1" ? 2 : 1;
        a.ps = ps;
        // 同一行的輸入欄位：E09 excl9 的 AD（bcftools TYPE=snp && N_ALT=1、GT 集合、FILTER 不含 PASS）與
        // E06 eps_vcf 的 AU/CU/GU/TU
        splitTabs(input, inF);
        if(inF.size() >= 10){
            std::vector<std::string> ifmt = splitChar(inF[8], ':');
            std::vector<std::string> ival = splitChar(inF[9], ':');
            auto get = [&](const char *key, std::string &out){
                auto it = std::find(ifmt.begin(), ifmt.end(), std::string(key));
                if(it == ifmt.end()){ return false; }
                size_t i = it - ifmt.begin();
                if(i >= ival.size()){ return false; }
                out = ival[i];
                return true;
            };
            std::string igt, iad;
            bool hasPass = false;
            for(const auto &f : splitChar(inF[6], ';')){ if(f == "PASS"){ hasPass = true; } }
            const bool gtOk = get("GT", igt) && (igt == "0/1" || igt == "1/0" || igt == "0|1" || igt == "1|0" || igt == "0|0");
            if(gtOk && !hasPass && get("AD", iad)){
                std::vector<std::string> ad = splitChar(iad, ',');
                long r, n2;
                if(ad.size() >= 2 && parseInt(ad[0], r) && parseInt(ad[1], n2)){
                    a.adKnown = true; a.adRef = r; a.adN = r + n2;
                }
            }
            long cc[4];
            bool ok = true;
            const char *keys[4] = {"AU", "CU", "GU", "TU"};
            for(int i = 0; i < 4 && ok; ++i){
                std::string v;
                ok = get(keys[i], v) && parseInt(v, cc[i]);
            }
            if(ok && a.ref >= 0 && a.alt >= 0){
                const long s = cc[0] + cc[1] + cc[2] + cc[3];
                if(s != 0){
                    a.baseKnown = true; a.baseTotal = s; a.baseThird = s - cc[a.ref] - cc[a.alt];
                }
            }
        }
        cs.anchors.push_back(a);
    });
}

// ---------------------------------------------------------------- BAM 計數（E03 count_hap.py，E09 版）
struct VoteAnchors {
    std::vector<int> pos; std::vector<int> ref, alt, hap; std::vector<long> ps;
};

// 已排序的 0-based ref 位置 → query 位置；非 M/=/X 為 -1（count_hap.qpos）
void queryPositions(const bam1_t *b, const int *refPos, int count, std::vector<int> &out){
    out.assign(count, -1);
    if(count == 0){ return; }
    const uint32_t *cigar = bam_get_cigar(b);
    const int nCigar = b->core.n_cigar;
    long rpos = b->core.pos;
    long qpos = 0;
    int i = 0;
    for(int c = 0; c < nCigar && i < count; ++c){
        const int op = bam_cigar_op(cigar[c]);
        const long len = bam_cigar_oplen(cigar[c]);
        const int type = bam_cigar_type(op);   // bit0: consumes query, bit1: consumes ref
        if(type & 2){
            const long rend = rpos + len;
            while(i < count && refPos[i] < rend){
                if(refPos[i] >= rpos && (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)){
                    out[i] = static_cast<int>(qpos + (refPos[i] - rpos));
                }
                ++i;
            }
            rpos = rend;
        }
        if(type & 1){ qpos += len; }
    }
}

int readBase(const bam1_t *b, int q){
    switch(bam_seqi(bam_get_seq(b), q)){ case 1: return 0; case 2: return 1; case 4: return 2; case 8: return 3; default: return -1; }
}

void countChromosome(const std::vector<std::string> &bamFiles, const std::string &fasta, const std::string &chrom,
                     const VoteAnchors &va, std::vector<Target> &targets, long &nReads){
    // block：最近錨點（排除同位置）的 PS，距離相同取左邊（count_hap）
    const std::vector<int> &apos = va.pos;
    for(Target &t : targets){
        t.tps = -1; t.dist = -1;
        if(apos.empty()){ continue; }
        const int p = t.pos0;
        const long k = std::lower_bound(apos.begin(), apos.end(), p) - apos.begin();
        long best = -1;
        for(long i : {k - 1, k, k + 1}){
            if(i < 0 || i >= static_cast<long>(apos.size()) || apos[i] == p){ continue; }
            if(best < 0){ best = i; continue; }
            const long di = std::labs(static_cast<long>(apos[i]) - p), db = std::labs(static_cast<long>(apos[best]) - p);
            if(di < db || (di == db && apos[i] < apos[best])){ best = i; }
        }
        if(best >= 0){ t.tps = va.ps[best]; t.dist = std::abs(apos[best] - p); }
    }
    std::vector<int> tpos(targets.size());
    for(size_t j = 0; j < targets.size(); ++j){ tpos[j] = targets[j].pos0; }

    std::vector<int> qa, qt;
    std::vector<int> vhap;
    std::vector<std::pair<long, std::pair<long, long>>> tally;
    for(const std::string &bamFile : bamFiles){
        samFile *fp = sam_open(bamFile.c_str(), "r");
        if(fp == nullptr){ std::cerr << "germline filter: cannot open " << bamFile << "\n"; exit(EXIT_FAILURE); }
        if(fp->format.format == cram){
            hts_set_fai_filename(fp, fasta.c_str());
            hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, 0x63E);
        }
        sam_hdr_t *hdr = sam_hdr_read(fp);
        hts_idx_t *idx = sam_index_load(fp, bamFile.c_str());
        const int tid = sam_hdr_name2tid(hdr, chrom.c_str());
        if(idx == nullptr || tid < 0){
            if(idx == nullptr){ std::cerr << "germline filter: cannot load index of " << bamFile << "\n"; exit(EXIT_FAILURE); }
            hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
            continue;
        }
        hts_itr_t *itr = sam_itr_queryi(idx, tid, 0, HTS_POS_MAX);
        bam1_t *b = bam_init1();
        while(sam_itr_next(fp, itr, b) >= 0){
            const int flag = b->core.flag;
            if((flag & kExcludeFlags) || (flag & kSupplementaryFlag) || b->core.qual < kMinMapq){ continue; }
            if(b->core.l_qseq == 0){ continue; }
            ++nReads;
            const long rs = b->core.pos;
            const long re = bam_endpos(b);
            const long ta = std::lower_bound(tpos.begin(), tpos.end(), rs) - tpos.begin();
            const long tb = std::lower_bound(tpos.begin(), tpos.end(), re) - tpos.begin();
            if(ta >= tb){ continue; }
            const long aa = std::lower_bound(apos.begin(), apos.end(), rs) - apos.begin();
            const long ab = std::lower_bound(apos.begin(), apos.end(), re) - apos.begin();
            const uint8_t *qual = bam_get_qual(b);
            // 錨點票：鹼基＝ALT → 帶 ALT 的 hap；＝REF → 另一條；其他或 BQ<10 不投
            queryPositions(b, apos.data() + aa, static_cast<int>(ab - aa), qa);
            vhap.assign(ab - aa, 0);
            tally.clear();
            for(long i = aa; i < ab; ++i){
                const int q = qa[i - aa];
                if(q < 0 || qual[q] < kMinBaseQuality){ continue; }
                const int base = readBase(b, q);
                int h = 0;
                if(base >= 0 && base == va.alt[i]){ h = va.hap[i]; }
                else if(base >= 0 && base == va.ref[i]){ h = 3 - va.hap[i]; }
                vhap[i - aa] = h;
                if(h == 0){ continue; }
                auto it = std::find_if(tally.begin(), tally.end(), [&](const auto &e){ return e.first == va.ps[i]; });
                if(it == tally.end()){ tally.push_back({va.ps[i], {0, 0}}); it = tally.end() - 1; }
                (h == 1 ? it->second.first : it->second.second) += 1;
            }
            queryPositions(b, tpos.data() + ta, static_cast<int>(tb - ta), qt);
            for(long j = ta; j < tb; ++j){
                Target &t = targets[j];
                const int q = qt[j - ta];
                int al;
                if(q < 0 || qual[q] < kMinBaseQuality){ al = 2; }
                else{
                    const int base = readBase(b, q);
                    al = (base >= 0 && base == t.refB) ? 0 : ((base >= 0 && base == t.altB) ? 1 : 2);
                }
                if(al == 2){ t.cnt[8] += 1; continue; }
                long v1 = 0, v2 = 0;
                if(t.tps >= 0){
                    for(const auto &e : tally){ if(e.first == t.tps){ v1 = e.second.first; v2 = e.second.second; break; } }
                    // 同位置的錨點不投票（vp != tpos）
                    auto lo = std::lower_bound(apos.begin() + aa, apos.begin() + ab, t.pos0);
                    for(auto it = lo; it != apos.begin() + ab && *it == t.pos0; ++it){
                        const long i = it - apos.begin();
                        if(va.ps[i] != t.tps){ continue; }
                        if(vhap[i - aa] == 1){ --v1; } else if(vhap[i - aa] == 2){ --v2; }
                    }
                }
                const long tot = v1 + v2, mx = std::max(v1, v2), mn = std::min(v1, v2);
                if(tot < 2){ t.cnt[4 + al] += 1; }
                else if(mx >= 2 && mx >= 2 * mn){
                    const int h = v1 > v2 ? 1 : 2;
                    t.cnt[(h == 1 ? 0 : 2) + al] += 1;
                    t.cnt[9] += tot; t.cnt[10] += 1;
                }
                else{ t.cnt[6 + al] += 1; }
            }
        }
        bam_destroy1(b);
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
    }
}

// ---------------------------------------------------------------- 逐 hap 摘要（decide6.hapcols）
struct HapSummary { int hapA; long k, n, kB, nB, kG, nG; bool noAlt; };
HapSummary summarize(const long *c){
    HapSummary s;
    const long k1 = c[1], k2 = c[3], n1 = c[0] + c[1], n2 = c[2] + c[3];
    const double fr1 = n1 > 0 ? static_cast<double>(k1) / n1 : -1;
    const double fr2 = n2 > 0 ? static_cast<double>(k2) / n2 : -1;
    s.hapA = k1 > k2 ? 1 : (k2 > k1 ? 2 : (fr2 > fr1 ? 2 : 1));
    s.k = s.hapA == 1 ? k1 : k2; s.n = s.hapA == 1 ? n1 : n2;
    s.kB = s.hapA == 1 ? k2 : k1; s.nB = s.hapA == 1 ? n2 : n1;
    s.noAlt = k1 == 0 && k2 == 0;
    s.kG = c[1] + c[3] + c[5] + c[7];
    s.nG = s.kG + c[0] + c[2] + c[4] + c[6];
    return s;
}

// ---------------------------------------------------------------- H（E03 decide.py，E06 只換 err）
std::vector<double> hCandidates(double p, double mj, double mn, double err){
    std::set<int> as;
    for(double x : {mj, mn}){
        if(std::isfinite(x) && roundHalfEven(x) >= 1){ as.insert(static_cast<int>(roundHalfEven(x))); }
    }
    if(as.empty()){ as.insert(1); }
    std::set<double> fs;
    for(int a : as){
        for(int m = 1; m <= a; ++m){ fs.insert((1 - err) * p * m / (p * a + 1 - p)); }
    }
    return std::vector<double>(fs.begin(), fs.end());
}
double hLr(long k, long n, double q, const std::vector<double> &fs){
    const double lg = binomLogpmf(k, n, q);
    double ls = kNegInf;
    for(double f : fs){ ls = std::max(ls, binomLogpmf(k, n, f)); }
    return lg - ls;
}
bool hJudgeable(long n, double q, const std::vector<double> &fs){
    // decide.judgeable：q 用 Python round(·,6)、fs 用 np.round(·,6)
    const double qr = pyRound6(q);
    std::vector<double> fr;
    for(double f : fs){ fr.push_back(npRound6(f)); }
    const double fmax = *std::max_element(fr.begin(), fr.end());
    if(n <= 0 || fmax >= qr){ return false; }
    double pKill = 0, pKeep = 0;
    for(long k = 0; k <= n; ++k){
        if(hLr(k, n, qr, fr) > 0){ pKill += binomPmf(k, n, fmax); }
        else{ pKeep += binomPmf(k, n, qr); }
    }
    return pKill <= kAlpha && pKeep <= kAlpha;
}

// ---------------------------------------------------------------- G（E06 decide6 g_states／g_lr／g_judge）
struct GStates { bool ok = false; std::vector<double> germ; double fmax = 0; };
GStates gStates(double p, double mj, double mn, double eps){
    GStates g;
    if(!(std::isfinite(p) && std::isfinite(mj) && std::isfinite(mn))){ return g; }
    int a = static_cast<int>(roundHalfEven(mj)), b = static_cast<int>(roundHalfEven(mn));
    const int hi = std::max(a, b), lo = std::min(a, b);
    a = hi; b = lo;
    const double d = p * (a + b) + 2 * (1 - p);
    if(d <= 0){ return g; }
    std::set<double> germ = {(p * a + 1 - p) / d, (p * b + 1 - p) / d, 1 - eps};
    for(double x : germ){ g.germ.push_back(npRound6(clip(x, eps, 1 - eps))); }
    const double fmax = p * std::max(a, 1) / d;
    g.fmax = pyRound6(std::min(std::max(fmax, eps), 1 - eps));
    g.ok = true;
    return g;
}
double gLr(long k, long n, const GStates &g){
    double lg = kNegInf;
    for(double x : g.germ){ lg = std::max(lg, binomLogpmf(k, n, x)); }
    const double fh = clip(std::min(static_cast<double>(k) / std::max(n, 1L), g.fmax), 1e-9, 1 - 1e-9);
    return lg - binomLogpmf(k, n, fh);
}
bool gJudgeable(long n, const GStates &g){
    double gLo = std::numeric_limits<double>::infinity();
    bool any = false;
    for(double x : g.germ){ if(x > g.fmax){ any = true; gLo = std::min(gLo, x); } }
    if(n <= 0 || !any){ return false; }
    double pKill = 0, pKeep = 0;
    for(long k = 0; k <= n; ++k){
        if(gLr(k, n, g) > 0){ pKill += binomPmf(k, n, g.fmax); }
        else{ pKeep += binomPmf(k, n, gLo); }
    }
    return pKill <= kAlpha && pKeep <= kAlpha;
}

// ---------------------------------------------------------------- gnomAD（E06 gnomad6.py）
class Gnomad {
    public:
        explicit Gnomad(const std::string &path){
            fp = hts_open(path.c_str(), "r");
            if(fp != nullptr){ tbx = tbx_index_load(path.c_str()); }
        }
        ~Gnomad(){
            if(tbx != nullptr){ tbx_destroy(tbx); }
            if(fp != nullptr){ hts_close(fp); }
            free(line.s);
        }
        bool ok() const { return fp != nullptr && tbx != nullptr; }
        // (chrom, pos1, REF, ALT) 相符紀錄的最大 AC；沒有命中回傳 -1
        long maxAc(const std::string &chrom, int pos1, const std::string &ref, const std::string &alt){
            const int tid = tbx_name2id(tbx, chrom.c_str());
            if(tid < 0){ return -1; }
            long best = -1;
            hts_itr_t *itr = tbx_itr_queryi(tbx, tid, pos1 - 1, pos1);
            if(itr == nullptr){ return -1; }
            std::vector<std::string> t;
            while(tbx_itr_next(fp, tbx, itr, &line) >= 0){
                splitTabs(std::string(line.s, line.l), t);
                if(t.size() < 8){ continue; }
                std::map<std::string, std::string> info;
                for(const auto &kv : splitChar(t[7], ';')){
                    size_t eq = kv.find('=');
                    if(eq != std::string::npos){ info[kv.substr(0, eq)] = kv.substr(eq + 1); }
                }
                auto acIt = info.find("AC");
                if(acIt == info.end()){ ++missingAc; continue; }
                std::vector<std::string> acs = splitChar(acIt->second, ',');
                std::vector<std::string> alts = splitChar(t[4], ',');
                const int recPos = std::atoi(t[1].c_str());
                for(size_t i = 0; i < alts.size(); ++i){
                    // 修剪共同後綴、再修剪前綴（gnomad6.norm）
                    std::string r = upper(t[3]), a = upper(alts[i]);
                    int q = recPos;
                    while(r.size() > 1 && a.size() > 1 && r.back() == a.back()){ r.pop_back(); a.pop_back(); }
                    while(r.size() > 1 && a.size() > 1 && r[0] == a[0]){ r.erase(0, 1); a.erase(0, 1); ++q; }
                    if(q == pos1 && r.size() == 1 && a.size() == 1 && r == ref && a == alt && i < acs.size()){
                        long ac;
                        if(parseInt(acs[i], ac)){ best = std::max(best, ac); }
                    }
                }
            }
            hts_itr_destroy(itr);
            return best;
        }
        long missingAc = 0;
    private:
        htsFile *fp = nullptr;
        tbx_t *tbx = nullptr;
        kstring_t line = {0, 0, nullptr};
};

void addRelabel(GermlineRelabelSet &relabel, const std::string &chrom, const Target &t){
    relabel.sites[chrom][t.pos0].push_back({t.ref, t.alt});
}

}  // namespace

// ---------------------------------------------------------------- 主流程
void run(const PhasingParameters &params, SnpParser &snpFile, ChrPhasingResult &chrPhasingResult,
         const Inputs &inputs, GermlineRelabelSet &relabel){
    const std::time_t begin = time(NULL);
    std::cerr << "germline filter: " << (inputs.highPurity
        ? "TDF > 0.9, H/G off, gnomAD AC>=2 on GT-selected SNVs"
        : "TDF <= 0.9, H+G with pseudo-het anchors excluded, gnomAD AC>=2 where both are undecidable") << "\n";
    if(!inputs.highPurity && (!inputs.purpleAvailable || inputs.copyNumbers == nullptr || inputs.copyNumbers->empty())){
        std::cerr << "germline filter: PURPLE result not available, filter not applied\n";
        return;
    }
    std::unique_ptr<Gnomad> gnomad;
    if(params.germlineGnomad.empty()){
        std::cerr << "germline filter: no gnomAD file (--germline-gnomad), gnomAD part skipped"
                  << (inputs.highPurity ? "; nothing to do at TDF > 0.9" : "; H/G only") << "\n";
        if(inputs.highPurity){ return; }
    }
    else{
        gnomad.reset(new Gnomad(params.germlineGnomad));
        if(!gnomad->ok()){
            std::cerr << "germline filter: cannot open gnomAD file or its tabix index: " << params.germlineGnomad << "\n";
            exit(EXIT_FAILURE);
        }
    }

    // ---- 輸出紀錄模擬
    std::map<std::string, ChromSites> sites;
    long nRecords = 0;
    collectSites(snpFile, params.snpFile, chrPhasingResult, inputs.highPurity, sites, nRecords);
    std::cerr << "germline filter: scanned " << nRecords << " records (" << difftime(time(NULL), begin) << "s)\n";

    std::ofstream dump;
    if(!params.germlineFilterDump.empty()){
        dump.open(params.germlineFilterDump);
        if(!dump){ std::cerr << "germline filter: cannot write " << params.germlineFilterDump << "\n"; exit(EXIT_FAILURE); }
    }

    // ---- TDF > 0.9：GT 挑到的 SNV 全部查 gnomAD
    if(inputs.highPurity){
        long n = 0, removed = 0;
        if(dump.is_open()){ dump << "chrom\tpos\tref\talt\tAC\tremoved\n"; }
        for(auto &entry : sites){
            for(const Target &t : entry.second.gtSelected){
                ++n;
                const long ac = gnomad->maxAc(entry.first, t.pos0 + 1, t.ref, t.alt);
                const bool rm = ac >= kMinGnomadAc;
                if(rm){ ++removed; addRelabel(relabel, entry.first, t); }
                if(dump.is_open()){ dump << entry.first << '\t' << t.pos0 + 1 << '\t' << t.ref << '\t' << t.alt << '\t' << ac << '\t' << rm << '\n'; }
            }
        }
        std::cerr << "germline filter: GT-selected SNVs " << n << ", gnomAD AC>=2 relabeled " << removed
                  << ", gnomAD records without AC " << gnomad->missingAc << " (" << difftime(time(NULL), begin) << "s)\n";
        return;
    }

    // ---- PURPLE：與 purity.tsv／cnv.somatic.tsv 相同的四位小數
    const double p = std::stod(purple::formatNumber(inputs.purplePurity));
    SegmentIndex segIndex;
    for(const auto &cn : *inputs.copyNumbers){
        segIndex.byChrom[cn.chromosome].push_back({cn.start, cn.end,
            std::stod(purple::formatNumber(purple::majorAlleleCopyNumber(cn))),
            std::stod(purple::formatNumber(purple::minorAlleleCopyNumber(cn)))});
    }
    for(auto &entry : segIndex.byChrom){
        std::stable_sort(entry.second.begin(), entry.second.end(), [](const Segment &a, const Segment &b){ return a.start < b.start; });
    }

    // ---- 偽 het 錨點（E09 excl9）、對照、投票錨點
    long nAnchors = 0, nAtPass = 0, nNoAd = 0, nFallback = 0, nExcluded = 0, nExcludedMain = 0;
    std::map<std::string, VoteAnchors> vote;
    for(auto &entry : sites){
        ChromSites &cs = entry.second;
        VoteAnchors &va = vote[entry.first];
        for(Anchor &a : cs.anchors){
            ++nAnchors;
            a.atPass = cs.passPos0.count(a.pos0) > 0;
            if(a.atPass){ ++nAtPass; continue; }
            if(!a.adKnown){ ++nNoAd; }
            else{
                double mjRaw, mnRaw;
                segIndex.find(entry.first, a.pos0 + 1, mjRaw, mnRaw);
                const int R = static_cast<int>(a.adRef), n = static_cast<int>(a.adN);
                const bool mainHom = binomTailGe(R, n, kEpsExclude) >= kAlphaHom;
                bool ok = std::isfinite(mjRaw) && std::isfinite(mnRaw);
                const double aa = std::max(0.0, roundHalfEven(ok ? mjRaw : 0.0));
                const double bb = std::max(0.0, roundHalfEven(ok ? mnRaw : 0.0));
                const double mj = std::max(aa, bb), mn = std::min(aa, bb);
                ok = ok && (mj + mn) > 0;
                bool hom;
                if(ok){
                    const double dd = p * (mj + mn) + 2 * (1 - p);
                    const double r1 = clip((p * mn + 1 - p) / dd, kEpsExclude, 1 - kEpsExclude);
                    const double r2 = clip((p * mj + 1 - p) / dd, kEpsExclude, 1 - kEpsExclude);
                    const double lh = binomLogpmf(R, n, kEpsExclude);
                    const double lhet = logHalfSumExp(binomLogpmf(R, n, r1), binomLogpmf(R, n, r2));
                    hom = (lh - lhet) > 0;
                }
                else{
                    ++nFallback;
                    hom = mainHom;
                }
                a.excluded = hom;
                if(hom){ ++nExcluded; }
                if(mainHom){ ++nExcludedMain; }
            }
            if(!a.excluded){
                va.pos.push_back(a.pos0); va.ref.push_back(a.ref); va.alt.push_back(a.alt);
                va.hap.push_back(a.hapAlt); va.ps.push_back(a.ps);
            }
        }
        // 對照：E03 的錨點（不排除 PASS 位置），chr1–22，段型非 NA。全部使用（E16 README §0.6）。
        for(const Anchor &a : cs.anchors){
            if(!isAutosome(entry.first)){ continue; }
            double mj, mn;
            segIndex.find(entry.first, a.pos0 + 1, mj, mn);
            if(segType(mj, mn) == SEG_NA){ continue; }
            Target t;
            t.pos0 = a.pos0;
            static const char *bases = "ACGT";
            t.ref = a.ref >= 0 ? std::string(1, bases[a.ref]) : "N";
            t.alt = a.alt >= 0 ? std::string(1, bases[a.alt]) : "N";
            t.refB = a.ref; t.altB = a.alt; t.isCtrl = true;
            cs.targets.push_back(t);
        }
        std::stable_sort(cs.targets.begin(), cs.targets.end(), [](const Target &x, const Target &y){ return x.pos0 < y.pos0; });
        // 投票錨點依位置排序（輸出 VCF 本來就排序；保險起見）
        std::vector<size_t> order(va.pos.size());
        for(size_t i = 0; i < order.size(); ++i){ order[i] = i; }
        std::stable_sort(order.begin(), order.end(), [&](size_t x, size_t y){ return va.pos[x] < va.pos[y]; });
        VoteAnchors sorted;
        for(size_t i : order){
            sorted.pos.push_back(va.pos[i]); sorted.ref.push_back(va.ref[i]); sorted.alt.push_back(va.alt[i]);
            sorted.hap.push_back(va.hap[i]); sorted.ps.push_back(va.ps[i]);
        }
        va = std::move(sorted);
    }
    std::cerr << "germline filter: anchors " << nAnchors << " (at PASS positions " << nAtPass << ", input AD unknown " << nNoAd
              << ", segment fallback " << nFallback << "), pseudo-het excluded " << nExcluded
              << " (E07 rule alone would exclude " << nExcludedMain << ")\n";

    // ---- BAM 計數
    std::vector<std::string> chroms;
    for(const auto &entry : sites){ if(!entry.second.targets.empty()){ chroms.push_back(entry.first); } }
    std::vector<long> readsPerChrom(chroms.size(), 0);
    const std::time_t countBegin = time(NULL);
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(int i = 0; i < static_cast<int>(chroms.size()); ++i){
        countChromosome(params.bamFile, params.fastaFile, chroms[i], vote.at(chroms[i]), sites.at(chroms[i]).targets, readsPerChrom[i]);
    }
    long nReads = 0;
    for(long r : readsPerChrom){ nReads += r; }
    std::cerr << "germline filter: BAM reads used " << nReads << " (" << difftime(time(NULL), countBegin) << "s)\n";

    // ---- ε 與 err（E06 v2：剔除疑似 hom-alt 對照）
    long baseTot = 0, baseThird = 0, baseN = 0;
    for(const auto &entry : sites){
        const ChromSites &cs = entry.second;
        if(!isAutosome(entry.first)){ continue; }
        for(const Anchor &a : cs.anchors){
            double mj, mn;
            segIndex.find(entry.first, a.pos0 + 1, mj, mn);
            if(segType(mj, mn) == SEG_NA || !a.baseKnown){ continue; }
            baseTot += a.baseTotal; baseThird += a.baseThird; ++baseN;
        }
    }
    const double epsRaw = static_cast<double>(baseThird) / (2.0 * std::max(baseTot, 1L));
    const double eps = std::max(epsRaw, kErrFloor);
    long wrong[4] = {0}, total[4] = {0}, used[4] = {0}, nCtrl[4] = {0}, nCtrlExcl[4] = {0};
    long wrong0[4] = {0}, total0[4] = {0};
    std::ofstream ctrlDump;
    if(dump.is_open()){
        ctrlDump.open(params.germlineFilterDump + ".ctrl.tsv");
        ctrlDump << "chrom\tpos\tref\talt\th1r\th1a\th2r\th2a\tunr\tuna\tcfr\tcfa\toth\tvote_sum\tn_voted\ttPS\tseg\tnG\tkG\thomalt_like\n";
    }
    for(const auto &entry : sites){
        for(const Target &t : entry.second.targets){
            if(!t.isCtrl){ continue; }
            double mj, mn;
            segIndex.find(entry.first, t.pos0 + 1, mj, mn);
            const int st = segType(mj, mn);
            const HapSummary s = summarize(t.cnt);
            const bool homLike = binomTailGe(static_cast<int>(s.nG - s.kG), static_cast<int>(s.nG), eps) >= kAlphaHom;
            const long wr = (s.n - s.k) + s.kB, tot = s.n + s.nB;
            for(int g : {st, 3}){
                ++nCtrl[g]; wrong0[g] += wr; total0[g] += tot;
                if(homLike){ ++nCtrlExcl[g]; continue; }
                wrong[g] += wr; total[g] += tot; ++used[g];
            }
            if(ctrlDump.is_open()){
                ctrlDump << entry.first << '\t' << t.pos0 + 1 << '\t' << t.ref << '\t' << t.alt;
                for(int c = 0; c < 11; ++c){ ctrlDump << '\t' << t.cnt[c]; }
                ctrlDump << '\t' << t.tps << '\t' << segName(st) << '\t' << s.nG << '\t' << s.kG << '\t' << homLike << '\n';
            }
        }
    }
    double err[4];
    for(int g = 0; g < 4; ++g){ err[g] = std::max(static_cast<double>(wrong[g]) / std::max(total[g], 1L), kErrFloor); }
    auto errOf = [&](int st){ return (st < 3 && used[st] >= kMinCtrlPerSeg) ? err[st] : err[3]; };
    std::cerr << "germline filter: eps " << eps << " (raw " << epsRaw << ", " << baseN << " controls)";
    for(int g = 0; g < 4; ++g){
        std::cerr << "; err " << (g == 3 ? "all" : segName(g)) << " " << err[g] << " (" << used[g] << "/" << nCtrl[g] << " controls)";
    }
    std::cerr << "\n";
    if(dump.is_open()){
        std::ofstream errDump(params.germlineFilterDump + ".err.tsv");
        errDump << "seg\tn_ctrl\tn_excl\terr_noexcl\terr_used\n";
        for(int g = 0; g < 4; ++g){
            errDump << (g == 3 ? "all" : segName(g)) << '\t' << nCtrl[g] << '\t' << nCtrlExcl[g] << '\t'
                    << std::max(static_cast<double>(wrong0[g]) / std::max(total0[g], 1L), kErrFloor) << '\t' << err[g] << '\n';
        }
        errDump << "eps\t" << baseN << "\t" << baseTot << "\t" << baseThird << "\t" << eps << '\n';
        std::ofstream exDump(params.germlineFilterDump + ".exclude.tsv");
        exDump << "chrom\tpos\n";
        for(const auto &entry : sites){
            for(const Anchor &a : entry.second.anchors){ if(a.excluded){ exDump << entry.first << '\t' << a.pos0 + 1 << '\n'; } }
        }
    }

    // ---- 判定（E09 decide9＋gnu9）
    long nPass = 0, nH = 0, nG = 0, nNone = 0, rmH = 0, rmG = 0, rmGn = 0, noneHit = 0;
    if(dump.is_open()){
        dump.precision(10);
        dump << "chrom\tpos\tref\talt\th1r\th1a\th2r\th2a\tunr\tuna\tcfr\tcfa\toth\tvote_sum\tn_voted\ttPS\tnear_anchor_dist"
                "\tMj\tMn\tseg\terr\tk\tn\tkG\tnG\tH_judge\tH_rm\tH_LR\tfmax_H\tG_judge\tG_rm\tG_LR\tfmax_G\tHG_path\tAC\tremoved\n";
    }
    for(const auto &entry : sites){
        for(const Target &t : entry.second.targets){
            if(t.isCtrl){ continue; }
            ++nPass;
            double mj, mn;
            segIndex.find(entry.first, t.pos0 + 1, mj, mn);
            const int st = segType(mj, mn);
            const double e = errOf(st);
            const double q = 1 - e;
            const HapSummary s = summarize(t.cnt);
            const std::vector<double> fs = hCandidates(p, mj, mn, e);
            const double fmaxH = *std::max_element(fs.begin(), fs.end());
            bool hj = false, hrm = false;
            double hlr = kNaN;
            if(s.n > 0 && !s.noAlt){
                hj = hJudgeable(s.n, q, fs);
                hlr = hLr(s.k, s.n, q, fs);
                hrm = hj && hlr > 0;
            }
            const GStates g = gStates(p, mj, mn, eps);
            bool gj = false, grm = false;
            double glr = kNaN;
            if(g.ok && s.nG > 0){
                gj = gJudgeable(s.nG, g);
                glr = gLr(s.kG, s.nG, g);
                grm = gj && glr > 0;
            }
            const char *path = hj ? "H" : (gj ? "G" : "none");
            bool rm = hj ? hrm : (gj ? grm : false);
            long ac = -1;
            if(hj){ ++nH; if(rm){ ++rmH; } }
            else if(gj){ ++nG; if(rm){ ++rmG; } }
            else{
                ++nNone;
                if(gnomad){
                    ac = gnomad->maxAc(entry.first, t.pos0 + 1, t.ref, t.alt);
                    if(ac >= 0){ ++noneHit; }
                    if(ac >= kMinGnomadAc){ rm = true; ++rmGn; }
                }
            }
            if(rm){ addRelabel(relabel, entry.first, t); }
            if(dump.is_open()){
                dump << entry.first << '\t' << t.pos0 + 1 << '\t' << t.ref << '\t' << t.alt;
                for(int c = 0; c < 11; ++c){ dump << '\t' << t.cnt[c]; }
                dump << '\t' << t.tps << '\t' << t.dist << '\t' << mj << '\t' << mn << '\t' << segName(st) << '\t' << e
                     << '\t' << s.k << '\t' << s.n << '\t' << s.kG << '\t' << s.nG
                     << '\t' << hj << '\t' << hrm << '\t' << hlr << '\t' << fmaxH
                     << '\t' << gj << '\t' << grm << '\t' << glr << '\t' << (g.ok ? g.fmax : kNaN)
                     << '\t' << path << '\t' << ac << '\t' << rm << '\n';
            }
        }
    }
    std::cerr << "germline filter: PASS SNVs " << nPass << "; H decidable " << nH << " (germline " << rmH
              << "), G decidable " << nG << " (germline " << rmG << "), neither " << nNone;
    if(gnomad){ std::cerr << " (gnomAD hit " << noneHit << ", AC>=2 " << rmGn << ")"; }
    std::cerr << "; relabeled " << (rmH + rmG + rmGn) << " as " << kGermlineFilterTag
              << " (" << difftime(time(NULL), begin) << "s)\n";
}

}

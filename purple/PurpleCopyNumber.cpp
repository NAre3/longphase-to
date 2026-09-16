#include "PurpleCopyNumber.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {
constexpr double EPS=1e-10;
bool zero(double x){return std::abs(x)<EPS;} bool le(double a,double b){return a-b<EPS;} bool gt(double a,double b){return a-b>EPS;}
double minor(const ObservedRegion&r){return r.tumorCopyNumber-r.tumorBaf*r.tumorCopyNumber;}
double major(const ObservedRegion&r){return r.tumorBaf*r.tumorCopyNumber;}

struct Combined {
    ObservedRegion value;
    CopyNumberMethod method=CopyNumberMethod::UNKNOWN;
    bool inferred=false,bafWeighted=true;
    int unweightedCount=1;
    std::vector<ObservedRegion> regions;
    Combined(bool weighted,const ObservedRegion&r):value(r),bafWeighted(weighted),regions{r}{if(r.germlineStatus!=GermlineStatus::DIPLOID)value.bafCount=0;}
    bool processed()const{return method!=CopyNumberMethod::UNKNOWN;}
    static double avg(long aw,double a,long bw,double b){if(zero(a))return b;return (aw*a+bw*b)/(aw+bw);}
    void extend(const ObservedRegion&r){
        const bool left=r.segment.start<=value.segment.start;
        value.segment.start=std::min(value.segment.start,r.segment.start); value.segment.end=std::max(value.segment.end,r.segment.end);
        value.segment.minStart=std::min(value.segment.minStart,r.segment.minStart); value.segment.maxStart=std::min(value.segment.maxStart,r.segment.maxStart);
        if(left){regions.insert(regions.begin(),r);value.segment.support=r.segment.support;value.segment.ratioSupport=r.segment.ratioSupport;}else regions.push_back(r);
    }
    void weightedExtend(const ObservedRegion&r){
        value.germlineStatus=GermlineStatus::DIPLOID;
        long aw=value.depthWindowCount,bw=r.depthWindowCount;
        if(!zero(r.tumorCopyNumber))value.tumorCopyNumber=avg(aw,value.tumorCopyNumber,bw,r.tumorCopyNumber);
        if(!zero(r.refNormalisedCopyNumber))value.refNormalisedCopyNumber=avg(aw,value.refNormalisedCopyNumber,bw,r.refNormalisedCopyNumber);
        if(!zero(r.gcContent))value.gcContent=avg(aw,value.gcContent,bw,r.gcContent);
        value.depthWindowCount+=r.depthWindowCount;
        aw=(bafWeighted&&(value.bafCount>0||r.bafCount>0))?value.bafCount:aw; bw=(bafWeighted&&(value.bafCount>0||r.bafCount>0))?r.bafCount:bw;
        if(!zero(r.observedBaf))value.observedBaf=avg(aw,value.observedBaf,bw,r.observedBaf);
        if(!zero(r.tumorBaf))value.tumorBaf=avg(aw,value.tumorBaf,bw,r.tumorBaf);
        value.bafCount+=r.bafCount; extend(r);
    }
    void inferBaf(double b){inferred=true;value.tumorBaf=b;value.bafCount=0;value.observedBaf=0;}
};

bool stop(const Combined&t,const ObservedRegion&n){return t.value.segment.start<n.segment.start?n.segment.support==SegmentSupport::CENTROMERE:t.value.segment.support==SegmentSupport::CENTROMERE;}
double relChange(double a,double b){double d=std::abs(std::max(a,b)-std::min(a,b));if(zero(d))return 0;if(zero(a)||zero(b))return 1;return d/std::abs(std::min(a,b));}
bool tolerance(const ObservedRegion&a,const ObservedRegion&b,double purity){
    double adjust=std::max(1.0,0.20/purity);int bc=std::min(a.bafCount,b.bafCount);
    if(bc>0){double maxCn=std::max(a.tumorCopyNumber,b.tumorCopyNumber);double tol=adjust*(0.3+0.5*maxCn/std::sqrt(bc));if(gt(std::abs(minor(a)-minor(b)),tol)&&gt(std::abs(a.observedBaf-b.observedBaf),0.03))return false;}
    int dc=std::min(a.depthWindowCount,b.depthWindowCount);double at=adjust*(0.3+1/std::sqrt(dc));double rt=adjust*(0.12+0.8/std::sqrt(dc));
    auto ok=[&](double x,double y){return le(std::abs(x-y),at)||le(relChange(x,y),rt);};return ok(a.tumorCopyNumber,b.tumorCopyNumber)||ok(a.refNormalisedCopyNumber,b.refNormalisedCopyNumber);
}
bool dubious(const ObservedRegion&r){return r.germlineStatus==GermlineStatus::DIPLOID&&r.depthWindowCount<150;}
bool valid(const ObservedRegion&r){return r.germlineStatus==GermlineStatus::DIPLOID&&r.depthWindowCount>=150;}
bool pushThrough(const std::vector<Combined>&v,int ti,int dir,double purity){int count=0;const Combined&t=v[ti];for(int i=ti+dir;i>=0&&i<(int)v.size();i+=dir){const auto&n=v[i].value;if(n.segment.start>t.value.segment.start&&n.segment.support==SegmentSupport::CENTROMERE)return count<150;bool d=dubious(n),in=tolerance(t.value,n,purity);if(d&&!in){count+=n.depthWindowCount;if(count>=150)return false;}if(valid(n))return in;if(n.segment.start<t.value.segment.start&&n.segment.support==SegmentSupport::CENTROMERE)return count<150;}return count<150;}
bool merge(std::vector<Combined>&v,int ti,int dir,double purity){Combined&t=v[ti];const ObservedRegion n=v[ti+dir].value;if(stop(t,n))return false;if(dubious(n)){if(tolerance(t.value,n,purity))t.weightedExtend(n);else if(pushThrough(v,ti,dir,purity))t.extend(n);else return false;}else if(!valid(n))t.extend(n);else if(tolerance(t.value,n,purity))t.weightedExtend(n);else return false;return true;}
int nextIndex(const std::vector<Combined>&v){int bi=-1,di=-1,bn=0,dn=0;for(int i=0;i<(int)v.size();++i)if(!v[i].processed()&&v[i].value.germlineStatus==GermlineStatus::DIPLOID){if(v[i].value.bafCount>bn){bn=v[i].value.bafCount;bi=i;}if(v[i].value.depthWindowCount>dn){dn=v[i].value.depthWindowCount;di=i;}}return bi>=0?bi:di;}
std::vector<Combined> extendDiploid(const std::vector<ObservedRegion>&input,double purity){bool bw=std::any_of(input.begin(),input.end(),[](auto&r){return r.bafCount>=50;});std::vector<Combined>v;for(auto&r:input)v.emplace_back(bw,r);for(int idx=nextIndex(v);idx>=0;idx=nextIndex(v)){v[idx].method=CopyNumberMethod::BAF_WEIGHTED;while(idx+1<(int)v.size()&&merge(v,idx,1,purity))v.erase(v.begin()+idx+1);while(idx>0&&merge(v,idx,-1,purity)){v.erase(v.begin()+idx-1);--idx;}}return v;}

bool acrocentric(const std::string&c){std::string x=lp::stripChrPrefix(c);return x=="13"||x=="14"||x=="15"||x=="21"||x=="22";}
void longArm(std::vector<Combined>&v){if(v.empty()||!acrocentric(v[0].value.segment.chromosome))return;int ci=-1;for(int i=0;i<(int)v.size();++i)if(v[i].value.segment.support==SegmentSupport::CENTROMERE){ci=i;break;}if(ci<=0||!v[ci].processed())return;double cn=v[ci].value.tumorCopyNumber;for(int i=ci-1;i>=0;){v[i].method=CopyNumberMethod::LONG_ARM;v[i].value.tumorCopyNumber=cn;if(v[i].value.segment.support!=SegmentSupport::NONE){--i;continue;}while(i>0){v[i].extend(v[i-1].value);v.erase(v.begin()+i-1);--i;if(v[i].value.segment.support!=SegmentSupport::NONE){--i;break;}}}}
void unknown(std::vector<Combined>&v){for(int i=0;i<(int)v.size();++i)if(!v[i].processed()){v[i].method=CopyNumberMethod::UNKNOWN;v[i].value.tumorCopyNumber=2;if(v[i].value.segment.support==SegmentSupport::NONE&&i>0&&v[i-1].method==CopyNumberMethod::UNKNOWN){v[i-1].extend(v[i].value);v.erase(v.begin()+i);--i;}}}

struct Gap{int ls=-1,lt=-1,rt=-1,rs=-1;bool valid()const{return ls>=0||rs>=0;}};
Gap nextGap(const std::vector<Combined>&v,bool cross){int ls=-1,lt=std::numeric_limits<int>::max();bool infer=false;for(int i=0;i<(int)v.size();++i){if(v[i].value.segment.support==SegmentSupport::CENTROMERE&&!cross){if(infer&&ls>-1)return{ls,lt,i-1,-1};ls=-1;infer=false;lt=std::numeric_limits<int>::max();}if(v[i].value.bafCount==0&&!v[i].inferred){infer=true;lt=std::min(lt,i);}else if(infer)return{ls,lt,i-1,i};else ls=i;}return infer?Gap{ls,lt,(int)v.size()-1,-1}:Gap{};}
double bafFor(double allele,double cn){if(le(cn,1))return 1;double ma=std::min(allele,cn),mi=cn-ma;return std::max(ma,mi)/(ma+mi);}
bool alleleDifferent(const ObservedRegion&a,const ObservedRegion&b,bool useMajor){return gt(std::abs((useMajor?major(a):minor(a))-(useMajor?major(b):minor(b))),0.5);}
bool eitherDifferent(const ObservedRegion&a,const ObservedRegion&b){return alleleDifferent(a,b,true)||alleleDifferent(a,b,false);}
const ObservedRegion* lookRight(const Gap&g,const std::vector<Combined>&v){const auto&src=v[g.rs].value;for(int i=g.rs+1;i<(int)v.size();++i){const auto&r=v[i].value;if(r.segment.support==SegmentSupport::CENTROMERE)return nullptr;if(r.bafCount>0&&r.segment.support!=SegmentSupport::CENTROMERE&&r.segment.support!=SegmentSupport::TELOMERE&&eitherDifferent(r,src))return &r;}return nullptr;}
const ObservedRegion* lookLeft(const Gap&g,const std::vector<Combined>&v){const auto&src=v[g.ls].value;for(int i=g.ls-1;i>=0;--i){if(v[i+1].value.segment.support==SegmentSupport::CENTROMERE)return nullptr;const auto&r=v[i].value;if(r.bafCount>0&&r.segment.support!=SegmentSupport::CENTROMERE&&r.segment.support!=SegmentSupport::TELOMERE&&eitherDifferent(r,src))return &r;}return nullptr;}
double closestAllele(double target,const ObservedRegion&s){return gt(std::abs(minor(s)-target),std::abs(major(s)-target))?major(s):minor(s);}
double movedAllele(const ObservedRegion&source,const ObservedRegion&primary,const ObservedRegion&secondary){bool md=alleleDifferent(primary,secondary,false),jd=alleleDifferent(primary,secondary,true);if(md&&jd){if(!gt(std::abs(minor(primary)-major(secondary)),0.5))return closestAllele((minor(primary)+major(secondary))/2,source);if(!gt(std::abs(major(primary)-minor(secondary)),0.5))return closestAllele((major(primary)+minor(secondary))/2,source);}return jd?minor(source):major(source);}
bool singleTinyGreater(const Gap&g,const std::vector<Combined>&v,const ObservedRegion&s){return g.lt==g.rt&&v[g.rt].value.segment.end-v[g.lt].value.segment.start+1<=30&&gt(v[g.lt].value.tumorCopyNumber,s.tumorCopyNumber);}
double targetAllele(const Gap&g,const std::vector<Combined>&v){
    if((g.ls>=0)^(g.rs>=0)){const ObservedRegion&s=g.ls>=0?v[g.ls].value:v[g.rs].value;if(singleTinyGreater(g,v,s))return minor(s);const ObservedRegion*d=g.ls>=0?lookLeft(g,v):lookRight(g,v);return d&&alleleDifferent(s,*d,false)?major(s):minor(s);}
    const ObservedRegion&left=v[g.ls].value,&right=v[g.rs].value;const ObservedRegion&primary=v[g.ls].value.bafCount>v[g.rs].value.bafCount?left:right;const ObservedRegion&secondary=&primary==&left?right:left;
    double minCn=std::numeric_limits<double>::max();for(int i=g.lt;i<=g.rt;++i)minCn=std::min(minCn,v[i].value.tumorCopyNumber);
    if(gt(major(left)-minCn,0.5)&&gt(major(right)-minCn,0.5))return minor(primary);
    if(eitherDifferent(primary,secondary))return movedAllele(primary,primary,secondary);
    if(singleTinyGreater(g,v,primary))return minor(primary);
    const ObservedRegion*rd=lookRight(g,v),*ld=lookLeft(g,v),*nearest=nullptr;
    if(!rd)nearest=ld;else if(!ld)nearest=rd;else{int dl=left.segment.end-ld->segment.start,dr=rd->segment.start-right.segment.start;nearest=dl<dr?ld:rd;}
    if(nearest){const ObservedRegion&other=nearest->segment.start>right.segment.start?right:left;return movedAllele(primary,*nearest,other);}
    return minor(primary);
}
void inferBafs(std::vector<Combined>&v){for(bool cross:{false,true})for(Gap g=nextGap(v,cross);g.valid();g=nextGap(v,cross)){double allele=targetAllele(g,v);for(int i=g.lt;i<=g.rt;++i)v[i].inferBaf(bafFor(allele,v[i].value.tumorCopyNumber));}}

std::vector<PurpleCopyNumber> convert(const std::vector<Combined>&v){std::vector<PurpleCopyNumber>o;for(int i=0;i<(int)v.size();++i){int start=v[i].value.segment.start;if(i>0){const auto&first=v[i].regions.front();if(first.germlineStatus==GermlineStatus::DIPLOID&&first.segment.minStart<first.segment.maxStart&&first.segment.minStart==o.back().end+1){start=(first.segment.minStart+first.segment.maxStart)/2;o.back().end=start-1;}}const auto&r=v[i].value;PurpleCopyNumber x;x.chromosome=r.segment.chromosome;x.start=start;x.end=r.segment.end;x.bafCount=r.bafCount;x.averageActualBaf=r.tumorBaf;x.averageObservedBaf=r.observedBaf;x.averageTumorCopyNumber=r.tumorCopyNumber;x.depthWindowCount=r.depthWindowCount;x.segmentStartSupport=r.segment.support;x.segmentEndSupport=i+1<(int)v.size()?v[i+1].value.segment.support:SegmentSupport::TELOMERE;x.method=v[i].method;x.gcContent=r.gcContent;x.minStart=r.segment.minStart;x.maxStart=r.segment.maxStart;o.push_back(x);}return o;}
}

std::vector<PurpleCopyNumber> buildCopyNumbers(const std::vector<ObservedRegion>&regions,const FittedPurity&fit){std::vector<PurpleCopyNumber>out;for(int chr=1;chr<=23;++chr){std::string key=chr==23?"X":std::to_string(chr);std::vector<ObservedRegion>cv;for(auto&r:regions)if(lp::stripChrPrefix(r.segment.chromosome)==key)cv.push_back(r);auto c=extendDiploid(cv,fit.purity);longArm(c);unknown(c);inferBafs(c);auto x=convert(c);out.insert(out.end(),x.begin(),x.end());}return out;}
void dumpCopyNumbers(const std::vector<PurpleCopyNumber>&v){std::vector<lp::CpDump::Row>rows;for(auto&x:v){double mi=le(x.averageActualBaf,0.5)&&std::abs(x.averageActualBaf-0.5)>EPS?0:std::max(0.0,(1-x.averageActualBaf)*x.averageTumorCopyNumber),ma=x.averageTumorCopyNumber-mi;rows.push_back({x.chromosome,x.start,x.chromosome+"\t"+std::to_string(x.start)+"\t"+std::to_string(x.end)+"\t"+std::to_string(x.bafCount)+"\t"+lp::CpDump::num(x.averageActualBaf)+"\t"+lp::CpDump::num(x.averageObservedBaf)+"\t"+lp::CpDump::num(x.averageTumorCopyNumber)+"\t"+std::to_string(x.depthWindowCount)+"\t"+segmentSupportName(x.segmentStartSupport)+"\t"+segmentSupportName(x.segmentEndSupport)+"\t"+copyNumberMethodName(x.method)+"\t"+lp::CpDump::num(x.gcContent)+"\t"+lp::CpDump::num(mi)+"\t"+lp::CpDump::num(ma)+"\t"+std::to_string(x.minStart)+"\t"+std::to_string(x.maxStart)+"\tfalse"});}lp::CpDump::write("CP-P8-copy-number","chromosome\tstart\tend\tbafCount\taverageActualBAF\taverageObservedBAF\taverageTumorCopyNumber\tdepthWindowCount\tsegmentStartSupport\tsegmentEndSupport\tmethod\tgcContent\tminorAlleleCopyNumber\tmajorAlleleCopyNumber\tminStart\tmaxStart\tsvSupport",rows);}
}

#include "SamRecordView.h"

#include "htslib/hts.h"

namespace amber {

SamRecordView::SamRecordView(const bam1_t *record)
    : mRecord(record)
{
    const bam1_core_t &core = record->core;

    mAlignmentStart = static_cast<int>(core.pos) + 1; // htslib 0-based → htsjdk 1-based
    mMappingQuality = core.qual;

    // 對應 SAMUtils.getAlignmentBlocks（SAMUtils.java:713-751），已與 jar bytecode 逐分支比對。
    // readBase 由 1 起算且 S 會前進、H 與 P 完全忽略 → read position 含 soft clip、不含 hard clip。
    const uint32_t *cigar = bam_get_cigar(record);
    int readBase = 1;
    int refBase = mAlignmentStart;
    int referenceLength = 0;

    mAlignmentBlocks.reserve(core.n_cigar);

    for(uint32_t i = 0; i < core.n_cigar; ++i){
        const int op = bam_cigar_op(cigar[i]);
        const int length = static_cast<int>(bam_cigar_oplen(cigar[i]));

        switch(op){
            case BAM_CHARD_CLIP: // H：忽略
            case BAM_CPAD:       // P：忽略
                break;
            case BAM_CSOFT_CLIP: // S：只消耗 query
            case BAM_CINS:       // I：只消耗 query
                readBase += length;
                break;
            case BAM_CREF_SKIP:  // N：只消耗 reference
            case BAM_CDEL:       // D：只消耗 reference
                refBase += length;
                referenceLength += length;
                break;
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                mAlignmentBlocks.push_back({readBase, refBase, length});
                readBase += length;
                refBase += length;
                referenceLength += length;
                break;
            default:
                // htsjdk 在此丟 IllegalStateException；BAM_CBACK 不會出現在合法的 BAM
                break;
        }
    }

    // htsjdk getAlignmentEnd()：unmapped 回傳 0，否則 start + referenceLength - 1。
    // referenceLength 為 0 時得到 start - 1，這正是與 bam_endpos() 分歧的那個情形。
    mAlignmentEnd = (core.flag & BAM_FUNMAP) ? 0 : mAlignmentStart + referenceLength - 1;
}

int SamRecordView::readPositionAtReferencePosition(int pos) const
{
    if(pos <= 0){
        return 0;
    }

    for(const AlignmentBlock &block : mAlignmentBlocks){
        // CoordMath.getEnd(start, length) = start + length - 1
        if(block.referenceStart + block.length - 1 >= pos){
            if(pos < block.referenceStart){
                // 中間隔了一段 deletion / ref skip。returnLastBaseIfDeleted 為 false → 0
                return 0;
            }
            return pos - block.referenceStart + block.readStart;
        }
        // returnLastBaseIfDeleted 為 false，故 htsjdk 的 lastAlignmentOffset 在此無作用
    }

    return 0; // 完全不與 read 重疊
}

int SamRecordView::referencePositionAtReadPosition(int position) const
{
    if(position == 0){
        return 0;
    }

    for(const AlignmentBlock &block : mAlignmentBlocks){
        if(block.readStart + block.length - 1 < position){
            continue;
        }
        if(position < block.readStart){
            return 0; // 落在 insertion 或 soft clip 之中
        }
        return block.referenceStart + position - block.readStart;
    }

    return 0;
}

char SamRecordView::baseAt(int readPosition) const
{
    const uint8_t *seq = bam_get_seq(mRecord);
    return seq_nt16_str[bam_seqi(seq, readPosition - 1)];
}

int SamRecordView::baseQualityAt(int readPosition) const
{
    const uint8_t *qual = bam_get_qual(mRecord);
    return qual[readPosition - 1];
}

}

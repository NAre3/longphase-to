#ifndef LP_SAMRECORDVIEW_H
#define LP_SAMRECORDVIEW_H

#include <string>
#include <vector>

#include "htslib/sam.h"

namespace lp {

// htsjdk SAMRecord 的座標轉換在 C++ 端的等價實作。
//
// 依據：RUN-004/behaviour-contract.md §3、§4（含 jar bytecode 交叉驗證）。
// 全部對外座標一律採 htsjdk 的約定：**1-based**，read position 含 soft clip、不含 hard clip，
// 亦即 read position 1 對應 bam_get_seq() 的索引 0。
//
// 每條 read 只建一次 alignment block 陣列（對應 htsjdk 在 SAMRecord 內快取 mAlignmentBlocks，
// SAMRecord.java:1786-1791）。PositionEvidenceChecker 對每個 (read, position) 都會呼叫
// getBaseQuality，若改成逐位點重走 CIGAR，熱路徑會被深度放大。
class SamRecordView
{
public:
    // 對應 htsjdk AlignmentBlock：三個欄位皆 1-based（AlignmentBlock.java:52-58）
    struct AlignmentBlock
    {
        int readStart;
        int referenceStart;
        int length;
    };

    explicit SamRecordView(const bam1_t *record);

    // 1-based inclusive，等同 htsjdk getAlignmentStart()
    int alignmentStart() const { return mAlignmentStart; }

    // 對應 htsjdk getAlignmentEnd()（SAMRecord.java:602-610）：
    //   alignmentStart + Cigar.getReferenceLength() - 1
    // 刻意**不使用 htslib 的 bam_endpos()**：bam_endpos 在 reference length 為 0 時
    // 強制 rlen = 1，會得到 alignmentStart，而 htsjdk 會得到 alignmentStart - 1。
    // 見 behaviour-contract.md §4.3。
    int alignmentEnd() const { return mAlignmentEnd; }

    int mappingQuality() const { return mMappingQuality; }

    // 對應 htsjdk getReadPositionAtReferencePosition(pos)，即
    // getReadPositionAtReferencePosition(rec, pos, returnLastBaseIfDeleted = false)
    // （SAMRecord.java:690-692 → 734-757；bytecode 已確認單參數多載傳入 false）。
    // 回傳 1-based read position；pos <= 0、pos 落在 D/N 之中、或 pos 在對齊範圍外皆回傳 0。
    int readPositionAtReferencePosition(int pos) const;

    // 對應 htsjdk getReferencePositionAtReadPosition(position)
    // （SAMRecord.java:643-646 → 665-678）。read position 落在 I 或 S 之中回傳 0。
    int referencePositionAtReadPosition(int position) const;

    // 1-based read position → 鹼基字元。對應 getReadString().charAt(readPosition - 1)。
    // 字元集為大寫，與 htsjdk 的 ACGTN= 相容（seq_nt16_str，hts.h:466）。
    char baseAt(int readPosition) const;

    // 1-based read position → binary phred。對應 getBaseQualities()[readPosition - 1]。
    // QUAL 與 SEQ 同索引基準（behaviour-contract.md §3.6）。
    int baseQualityAt(int readPosition) const;

private:
    const bam1_t *mRecord;
    int mAlignmentStart;
    int mAlignmentEnd;
    int mMappingQuality;
    std::vector<AlignmentBlock> mAlignmentBlocks;
};

}

#endif

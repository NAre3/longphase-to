#ifndef LP_CPDUMP_H
#define LP_CPDUMP_H

#include <string>
#include <vector>

namespace lp {

// 保真度驗證用的 checkpoint dump（對應 Java 端插入的同名 dump）。格式：
// TSV、\t 分隔、\n 行尾、固定欄名；per-locus 由 dump 端依 (chromosome, position) 穩定排序。
// 比對端以數值解析後比較，因此浮點的文字形式不需與 Java 逐字相同。
class CpDump
{
public:
    struct Row
    {
        std::string chromosome;
        int position;
        std::string line;
    };

    static bool enabled();
    static void setDir(const std::string &dir);

    // 依 (chromosome, position) 穩定排序後輸出
    static void writeSorted(const std::string &checkpoint, const std::string &header, std::vector<Row> &rows);

    // 保留呼叫端順序輸出
    static void write(const std::string &checkpoint, const std::string &header, const std::vector<Row> &rows);

    static std::string num(double value);
};

}

#endif

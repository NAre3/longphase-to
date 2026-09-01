#ifndef AMBER_AMBEROUTPUT_H
#define AMBER_AMBEROUTPUT_H

#include <string>
#include <vector>

namespace amber {

// 對應 hmf-common 的 AmberBAF（AmberBAF.java）。
struct AmberBAF
{
    std::string chromosome;
    int position = 0;
    double tumorBAF = 0.0;
    int tumorDepth = 0;
    double normalBAF = 0.0;
    int normalDepth = 0;

    // AmberBAF.java:tumorModifiedBAF / normalModifiedBAF
    double tumorModifiedBAF() const;
    double normalModifiedBAF() const;
};

// 對應 AmberBAFFile.write（4 位小數、非有限值寫 "0"）。以 gzip 寫出。
void writeAmberBafFile(const std::string &path, const std::vector<AmberBAF> &bafs);

// 對應 AmberQCFile.write + AmberQC.status()。
void writeAmberQcFile(const std::string &path, double contamination, double consanguinityProportion);

}

#endif

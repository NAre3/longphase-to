#include "HumanChromosome.h"

#include <cctype>

namespace amber {

std::string stripChrPrefix(const std::string &chromosome){
    if(chromosome.size() >= 3){
        if(std::tolower(chromosome[0]) == 'c' && std::tolower(chromosome[1]) == 'h' && std::tolower(chromosome[2]) == 'r'){
            return chromosome.substr(3);
        }
    }
    return chromosome;
}

bool isHumanChromosome(const std::string &chromosome){
    const std::string trimmed = stripChrPrefix(chromosome);

    if(trimmed.empty()){
        return false;
    }

    // Java 端用 StringUtils.isNumeric：全部字元皆為數字才算數值
    bool numeric = true;
    for(char c : trimmed){
        if(!std::isdigit(static_cast<unsigned char>(c))){
            numeric = false;
            break;
        }
    }

    if(numeric){
        // 位數過多時 stoi 會丟例外，先擋掉；Java 端的來源檔不會出現這種值
        if(trimmed.size() > 9){
            return false;
        }
        const int value = std::stoi(trimmed);
        return value >= 1 && value <= 22;
    }

    return trimmed == "X" || trimmed == "Y";
}

}

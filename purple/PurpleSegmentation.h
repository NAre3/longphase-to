#ifndef PURPLE_SEGMENTATION_H
#define PURPLE_SEGMENTATION_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

std::vector<SupportSegment> createSupportSegments(const InputData &data, const std::string &referenceFasta);
void dumpSupportSegments(const std::vector<SupportSegment> &segments);

}

#endif

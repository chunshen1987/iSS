#ifndef SRC_URQMD_OUTPUT_H_
#define SRC_URQMD_OUTPUT_H_

#include <string>
#include <vector>

#include "data_struct.h"

namespace iSS_data {

void write_urqmd_samples(const std::string &filename,
                         const std::string &oscar_header_filename,
                         const std::vector<std::vector<iSS_Hadron>*> &events,
                         std::size_t first_event = 0);

}  // namespace iSS_data

#endif  // SRC_URQMD_OUTPUT_H_

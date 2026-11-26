#ifndef SRC_SAMPLE_OUTPUT_H_
#define SRC_SAMPLE_OUTPUT_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "data_struct.h"
#include "zlib.h"

namespace iSS_data {

[[noreturn]] inline void sample_output_error(const std::string &filename,
                                             const char *reason) {
    std::cerr << "[Error] Sample output " << filename << ": " << reason << '\n';
    std::exit(1);
}

inline void write_binary_samples(
        const std::string &filename,
        const std::vector<std::vector<iSS_Hadron>*> &events) {
    std::ofstream output(filename, std::ios::binary);
    if (!output) sample_output_error(filename, "cannot open file");
    for (const auto *event : events) {
        const int count = event->size();
        output.write(reinterpret_cast<const char*>(&count), sizeof(count));
        for (const auto &particle : *event) {
            const float values[] = {
                particle.mass, particle.t, particle.x, particle.y, particle.z,
                particle.E, particle.px, particle.py, particle.pz,
            };
            output.write(reinterpret_cast<const char*>(&particle.pid), sizeof(int));
            output.write(reinterpret_cast<const char*>(values), sizeof(values));
        }
        if (!output) sample_output_error(filename, "write failed");
    }
    output.close();
    if (!output) sample_output_error(filename, "write failed");
}

inline void write_gzip_samples(
        const std::string &filename,
        const std::vector<std::vector<iSS_Hadron>*> &events) {
    gzFile output = gzopen(filename.c_str(), "wb");
    if (!output) sample_output_error(filename, "cannot open file");
    for (const auto *event : events) {
        if (gzprintf(output, "%d \n", static_cast<int>(event->size())) <= 0) {
            gzclose(output);
            sample_output_error(filename, "write failed");
        }
        for (const auto &particle : *event) {
            if (gzprintf(output, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
                         particle.pid, particle.mass, particle.t, particle.x,
                         particle.y, particle.z, particle.E, particle.px,
                         particle.py, particle.pz) <= 0) {
                gzclose(output);
                sample_output_error(filename, "write failed");
            }
        }
    }
    if (gzclose(output) != Z_OK) sample_output_error(filename, "write failed");
}

}  // namespace iSS_data

#endif  // SRC_SAMPLE_OUTPUT_H_

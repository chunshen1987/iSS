#ifndef SRC_SAMPLE_OUTPUT_H_
#define SRC_SAMPLE_OUTPUT_H_

#include <cstdlib>
#include <cstring>
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

class BinarySampleWriter {
 public:
    explicit BinarySampleWriter(const std::string &filename)
        : filename_(filename), output_(filename, std::ios::binary),
          buffer_(4096*record_bytes) {
        if (!output_) sample_output_error(filename_, "cannot open file");
    }

    void write(const std::vector<iSS_Hadron> &event) {
        const int count = event.size();
        output_.write(reinterpret_cast<const char*>(&count), sizeof(count));
        std::size_t buffered = 0;
        for (const auto &particle : event) {
            const float values[] = {
                particle.mass, particle.t, particle.x, particle.y, particle.z,
                particle.E, particle.px, particle.py, particle.pz,
            };
            char *record = buffer_.data() + buffered;
            std::memcpy(record, &particle.pid, sizeof(int));
            std::memcpy(record + sizeof(int), values, sizeof(values));
            buffered += record_bytes;
            if (buffered == buffer_.size()) {
                output_.write(buffer_.data(), buffered);
                buffered = 0;
            }
        }
        if (buffered != 0) output_.write(buffer_.data(), buffered);
        if (!output_) sample_output_error(filename_, "write failed");
    }

    void flush() {
        if (!output_.flush()) sample_output_error(filename_, "write failed");
    }

    void close() {
        output_.close();
        if (!output_) sample_output_error(filename_, "write failed");
    }

 private:
    // Serialize explicitly, independent of iSS_Hadron padding.
    static constexpr std::size_t record_bytes = sizeof(int) + 9*sizeof(float);
    std::string filename_;
    std::ofstream output_;
    std::vector<char> buffer_;
};

class GzipSampleWriter {
 public:
    explicit GzipSampleWriter(const std::string &filename)
        : filename_(filename), output_(gzopen(filename.c_str(), "wb")) {
        if (!output_) sample_output_error(filename_, "cannot open file");
    }

    ~GzipSampleWriter() {
        if (output_) gzclose(output_);
    }

    GzipSampleWriter(const GzipSampleWriter &) = delete;
    GzipSampleWriter &operator=(const GzipSampleWriter &) = delete;

    void write(const std::vector<iSS_Hadron> &event) {
        if (gzprintf(output_, "%d \n", static_cast<int>(event.size())) <= 0) fail();
        for (const auto &particle : event) {
            if (gzprintf(output_, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
                         particle.pid, particle.mass, particle.t, particle.x,
                         particle.y, particle.z, particle.E, particle.px,
                         particle.py, particle.pz) <= 0) fail();
        }
    }

    void flush() {
        if (gzflush(output_, Z_SYNC_FLUSH) != Z_OK) fail();
    }

    void close() {
        const int result = gzclose(output_);
        output_ = nullptr;
        if (result != Z_OK) sample_output_error(filename_, "write failed");
    }

 private:
    [[noreturn]] void fail() {
        gzclose(output_);
        sample_output_error(filename_, "write failed");
    }

    std::string filename_;
    gzFile output_;
};

inline void write_binary_samples(
        const std::string &filename,
        const std::vector<std::vector<iSS_Hadron>*> &events) {
    BinarySampleWriter output(filename);
    for (const auto *event : events) output.write(*event);
    output.close();
}

inline void write_gzip_samples(
        const std::string &filename,
        const std::vector<std::vector<iSS_Hadron>*> &events) {
    GzipSampleWriter output(filename);
    for (const auto *event : events) output.write(*event);
    output.close();
}

}  // namespace iSS_data

#endif  // SRC_SAMPLE_OUTPUT_H_

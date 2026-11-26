#include <iostream>
#include <vector>

#include "iSS.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: TestRepeatedSampling <MUSIC work directory>\n";
        return 1;
    }
    iSS sampler(argv[1], "iSS_tables", "iSS_tables", "iSS_parameters.dat",
                "surface.dat");
    auto *parameters = sampler.paraRdr_ptr;
    parameters->setVal("regulateEOS", 0);
    parameters->setVal("include_deltaf_shear", 0);
    parameters->setVal("include_deltaf_bulk", 0);
    parameters->setVal("use_binary_format", 0);
    parameters->setVal("number_of_repeated_sampling", 1);
    parameters->setVal("randomSeed", 1);

    if (sampler.shell() != 0) return 1;
    const auto first_event = *sampler.get_hadron_list_iev(0);
    if (sampler.shell() != 0) return 1;
    const auto &second_event = *sampler.get_hadron_list_iev(0);
    if (first_event.empty() || first_event.size() != second_event.size()) {
        std::cerr << "Repeating shell() changed the sampled multiplicity\n";
        return 1;
    }
    for (std::size_t i = 0; i < first_event.size(); ++i) {
        const auto &first = first_event[i];
        const auto &second = second_event[i];
        if (first.pid != second.pid || first.E != second.E || first.px != second.px
                || first.py != second.py || first.pz != second.pz) {
            std::cerr << "Repeating shell() changed the sampled momenta\n";
            return 1;
        }
    }
    std::cout << "Repeated shell() calls produce the same fixed-seed event.\n";
}

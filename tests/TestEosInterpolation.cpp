#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "eos_4D.h"

namespace {

// Header order is B, Q, S, T. The payload varies fastest along B.
const std::array<float, 12> grid = {{
    -0.5f, -0.25f, -0.75f, 0.125f,
     0.5f,  0.5f,   0.5f,  0.125f,
     3.f,   2.f,    2.f,   3.f,
}};

// A separately linear function in all four axes is exactly representable by
// multilinear interpolation. Distinct slopes and a cross term expose axis
// swaps, wrong strides, and missing corner weights.
double table_value(int field, double t, double b, double q, double s) {
    return 12. + field + 2.*t + 3.*b + 5.*q + 7.*s + 8.*t*b*q*s;
}

void write_table(const std::string &path, int first_field, int field_count,
                 bool has_header) {
    std::vector<float> values;
    if (has_header) values.assign(grid.begin(), grid.end());
    for (int t = 0; t < grid[11]; ++t) {
        for (int s = 0; s < grid[10]; ++s) {
            for (int q = 0; q < grid[9]; ++q) {
                for (int b = 0; b < grid[8]; ++b) {
                    for (int field = first_field;
                         field < first_field + field_count; ++field) {
                        values.push_back(table_value(
                            field, grid[3] + t*grid[7], grid[0] + b*grid[4],
                            grid[1] + q*grid[5], grid[2] + s*grid[6]));
                    }
                }
            }
        }
    }

    INFO(path);
    std::ofstream output(path, std::ios::binary);
    REQUIRE(output.is_open());
    output.write(reinterpret_cast<const char *>(values.data()),
                 values.size()*sizeof(float));
    output.close();
    REQUIRE(output.good());
}

std::array<double, 4> physical_input(const std::array<double, 4> &coordinates) {
    // Invert the EoS table coordinate map to obtain (e, nB, nQ, nS).
    const double t = coordinates[0];
    const double b = coordinates[1];
    const double q = coordinates[2];
    const double s = coordinates[3];
    const double hbarc_cubed = 0.19733*0.19733*0.19733;
    const double density_scale = t*t/(3.*hbarc_cubed);
    const double pi = std::acos(-1.);
    return {{19.*pi*pi*t*t*t*t/(12.*hbarc_cubed),
             (b - s)*density_scale,
             (2.*q + s)*density_scale,
             (3.*s + q - b)*density_scale}};
}

}  // namespace

TEST_CASE("4D EoS lookups reproduce a multilinear table and clamp its edges") {
    const std::string eos_dir = "iSS_tables/EOS_tables/HRG4D/";
    const std::string deltaf_dir = "iSS_tables/deltaf_tables/urqmd/";
    const char *fields[] = {"p", "t", "mub", "mus", "muq"};
    for (int i = 0; i < 5; ++i) {
        write_table(eos_dir + "HRG_" + fields[i] + "_b.dat", i, 1, true);
    }
    write_table(deltaf_dir + "NEoS4D_22mom_deltafCoeff.bin", 5, 6, false);
    write_table(deltaf_dir + "NEoS4D_CE_deltafCoeff.bin", 11, 3, false);

    const struct {
        std::array<double, 4> input;
        std::array<double, 4> expected;
    } cases[] = {
        {{0.31, 0.20, -0.07, -0.40}, {0.31, 0.20, -0.07, -0.40}},
        {{0.25, 0.00,  0.25, -0.75}, {0.25, 0.00,  0.25, -0.75}},
        {{0.50, 0.70, -0.40, -1.00}, {0.375, 0.50, -0.25, -0.75}},
        {{0.08,-0.70,  0.40,  0.00}, {0.125,-0.50,  0.25, -0.25}},
    };
    EOS_4D eos;
    eos.initialize_eos();
    for (int deltaf_type = -1; deltaf_type <= 1; ++deltaf_type) {
        CAPTURE(deltaf_type);
        if (deltaf_type >= 0) eos.initialize_dfCoeffs(deltaf_type);
        const std::size_t coefficient_counts[] = {0, 6, 3};
        const std::size_t coefficient_count = coefficient_counts[deltaf_type + 1];
        const int first_coefficient = deltaf_type == 1 ? 11 : 5;

        for (const auto &test : cases) {
            const auto input = physical_input(test.input);
            INFO("Grid coordinates T, B, Q, S: " << test.input[0] << ", "
                 << test.input[1] << ", " << test.input[2] << ", "
                 << test.input[3]);
            std::vector<float> thermal;
            std::vector<double> deltaf;
            std::vector<float> combined_thermal;
            std::vector<double> combined_deltaf;
            eos.getThermalVariables(input[0], input[1], input[2], input[3], thermal);
            eos.getDeltafCoeffs(input[0], input[1], input[2], input[3], deltaf);
            eos.getThermalVariablesAndDeltafCoeffs(
                input[0], input[1], input[2], input[3], combined_thermal,
                combined_deltaf);
            REQUIRE(thermal.size() == 5);
            REQUIRE(combined_thermal.size() == thermal.size());
            REQUIRE(deltaf.size() == coefficient_count);
            REQUIRE(combined_deltaf.size() == coefficient_count);

            // The tables and weights use float. Allow rounding in a sum of
            // 16 weighted corners as well as the round trip through (e, n).
            const double tolerance = 2.e-5;
            const auto &point = test.expected;
            for (std::size_t i = 0; i < thermal.size(); ++i) {
                CAPTURE(i);
                const double expected = table_value(i, point[0], point[1],
                                                    point[2], point[3]);
                CHECK(std::abs(thermal[i] - expected) < tolerance);
                CHECK(std::abs(combined_thermal[i] - expected) < tolerance);
            }
            for (std::size_t i = 0; i < coefficient_count; ++i) {
                CAPTURE(i);
                const double expected = table_value(
                    first_coefficient + i, point[0], point[1], point[2], point[3]);
                CHECK(std::abs(deltaf[i] - expected) < tolerance);
                CHECK(std::abs(combined_deltaf[i] - expected) < tolerance);
            }
        }
    }
}

#include <memory>
#include <vector>

#include "RandomVariable1DArray.h"
#include "RandomVariable2DArray.h"
#include "doctest.h"

TEST_CASE("One-dimensional sampling selects the correct CDF interval") {
    auto rng = std::make_shared<RandomUtil::Random>(17L, 0., 1.);
    std::vector<double> weights{2., 0., 3., 5.};
    RandomVariable1DArray sampler(&weights, rng, 0.);
    CHECK(sampler.return_sum() == 10.);

    // Keep the draw entirely within one interval, with probes on both sides
    // of each boundary. The zero-weight bin has no interval of its own.
    RandomUtil::Random below_first_boundary(17L, 0.199, 0.1999);
    RandomUtil::Random above_first_boundary(17L, 0.2001, 0.201);
    RandomUtil::Random below_second_boundary(17L, 0.499, 0.4999);
    RandomUtil::Random above_second_boundary(17L, 0.5001, 0.501);
    RandomUtil::Random near_one(17L, 0.999, 1.);
    CHECK(sampler.rand(below_first_boundary) == 0);
    CHECK(sampler.rand(above_first_boundary) == 2);
    CHECK(sampler.rand(below_second_boundary) == 2);
    CHECK(sampler.rand(above_second_boundary) == 3);
    CHECK(sampler.rand(near_one) == 3);
}

TEST_CASE("Two-dimensional sampling preserves row and column order") {
    double first_row[] = {1., 0., 2.};
    double second_row[] = {3., 0., 4.};
    double *weights[] = {first_row, second_row};
    const struct {
        double draw;
        long row;
        long column;
    } cases[] = {{0.05, 0, 0}, {0.2, 0, 2}, {0.4, 1, 0}, {0.8, 1, 2}};

    for (const auto &test : cases) {
        CAPTURE(test.draw);
        auto rng = std::make_shared<RandomUtil::Random>(17L, test.draw,
                                                        test.draw + 0.01);
        RandomVariable2DArray sampler(weights, 3, 2, rng, 0.);
        long row = -1;
        long column = -1;
        sampler.sampleAccToInvCDF(&row, &column);
        CHECK(row == test.row);
        CHECK(column == test.column);
    }
}

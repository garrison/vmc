#include "Lattice.hpp"

#include <gtest/gtest.h>

TEST(Lattice, Basic) {
    lw_vector<int, MAX_DIMENSION> dimensions;
    dimensions.push_back(8);
    dimensions.push_back(4);
    Lattice lattice(dimensions);

    EXPECT_EQ(lattice.n_dimensions(), 2);
    EXPECT_EQ(lattice.basis_indices, 1);
    EXPECT_EQ(lattice.dimensions[0], 8);
    EXPECT_EQ(lattice.dimensions[1], 4);
}

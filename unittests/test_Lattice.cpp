#include "Lattice.hpp"

#include <gtest/gtest.h>

TEST(Lattice, VectorConstructor) {
    Lattice::DimensionVector dimensions;
    dimensions.push_back(8);
    dimensions.push_back(4);
    Lattice lattice(dimensions);

    EXPECT_EQ(lattice.n_dimensions(), 2);
    EXPECT_EQ(lattice.basis_indices, 1);
    EXPECT_EQ(lattice.dimensions[0], 8);
    EXPECT_EQ(lattice.dimensions[1], 4);
}

TEST(Lattice, InitializerListConstructor) {
    Lattice lattice({ 8, 4 });

    EXPECT_EQ(lattice.n_dimensions(), 2);
    EXPECT_EQ(lattice.basis_indices, 1);
    EXPECT_EQ(lattice.dimensions[0], 8);
    EXPECT_EQ(lattice.dimensions[1], 4);
}

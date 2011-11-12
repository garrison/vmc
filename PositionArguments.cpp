#include "PositionArguments.hpp"

PositionArguments::PositionArguments (const std::vector<unsigned int> &r_, unsigned int N_sites)
    : r(r_),
      positions(N_sites)
{
    BOOST_ASSERT(get_N_filled() <= get_N_sites());
    for (unsigned int i = 0; i < r.size(); ++i) {
	BOOST_ASSERT(r[i] < get_N_sites()); // ensure valid site
	++positions[r[i]];
	BOOST_ASSERT(positions[r[i]] == 1); // no double filling allowed
    }
}

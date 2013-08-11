#include "PositionArguments.hpp"

PositionArguments::PositionArguments (const std::vector<std::vector<unsigned int> > &r_, unsigned int N_sites)
    : r(r_),
      positions(r.size())
{
    assert(r_.size() > 0);

    _populate_positions(N_sites);
}

void PositionArguments::reset (const std::vector<std::vector<unsigned int> > &r_)
{
    assert(r_.size() > 0);

#ifndef NDEBUG
    // this condition is not strictly necessary, of course, but for now it
    // helps assure we are using things correctly.  feel free to remove later.
    assert(r.size() == r_.size());
    for (unsigned int i = 0; i < get_N_species(); ++i)
        assert(r[i].size() == r_[i].size());
#endif

    r = r_;
    _populate_positions(get_N_sites());
}

int PositionArguments::particle_index_at_position (unsigned int position, unsigned int species) const
{
    assert(position < get_N_sites());
    assert(species < get_N_species());

    if (positions[species][position] == 0)
        return -1;
    for (unsigned int i = 0; i < get_N_filled(species); ++i) {
        if (r[species][i] == position)
            return i;
    }

    // should not be reached
    assert(false);
    return -1;
}

void PositionArguments::_populate_positions (unsigned int N_sites)
{
    N_filled_total = 0;

    // reset positions vectors
    for (unsigned int i = 0; i < get_N_species(); ++i)
        positions[i].assign(N_sites, 0);

    // populate them, while also determining N_filled_total
    for (unsigned int species = 0; species < r.size(); ++species) {
        // we don't allow double occupancy
        assert(get_N_filled(species) < N_sites);
        for (unsigned int j = 0; j < r[species].size() ; ++j) {
            assert(r[species][j] < N_sites);
            ++positions[species][r[species][j]];
            assert(positions[species][r[species][j]] == 1);
        }
        N_filled_total += get_N_filled(species);
    }
}

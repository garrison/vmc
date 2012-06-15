#include <set>
#include <utility>

#include "Move.hpp"

bool Move::is_valid_for (const PositionArguments &r) const
{
    const unsigned int N_sites = r.get_N_sites();

    // the pair represents (species, site_index).  This should be unique among
    // all things involved in the move.
    std::set<std::pair<unsigned int, unsigned int> > unique_pair_set;
    for (Move::const_iterator i = this->begin(); i != this->end(); ++i) {
        if (!(r.particle_is_valid(i->particle)
              && i->destination < N_sites
              && r[i->particle] != i->destination // null single particle moves are not allowed.
              && !r.is_occupied(i->destination, i->particle.species) // enforce PAULI
              && unique_pair_set.insert(std::make_pair(i->particle.species, r[i->particle])).second
              && unique_pair_set.insert(std::make_pair(i->particle.species, i->destination)).second))
            return false;
    }
    return true;
}

#ifndef _POSITION_ARGUMENTS_HPP
#define _POSITION_ARGUMENTS_HPP

#include <vector>

#include <boost/assert.hpp>

/**
 * Represents a particle in a system
 */
struct Particle
{
    unsigned int index;
    unsigned int species;

    Particle (void)
        {
        }

    Particle (unsigned int index_, unsigned int species_)
        : index(index_),
          species(species_)
        {
        }
};

/**
 * Stores each particle's position such that we can also quickly query to see
 * whether a given site is vacant or not.
 *
 * The assertions only allow one particle per site, but this restriction could
 * easily be removed.
 */
class PositionArguments
{
public:
    /**
     * Constructor
     *
     * @param r_ a vector of a vector which, for each species, contains the
     * positions of each particle
     */
    PositionArguments (const std::vector<std::vector<unsigned int> > &r_, unsigned int N_sites);

    /**
     * Resets the positions according to the given vector
     */
    void reset (const std::vector<std::vector<unsigned int> > &r_);

    /**
     * Returns the site index of the given particle
     */
    unsigned int operator[] (Particle particle) const
        {
            BOOST_ASSERT(particle_is_valid(particle));
            return r[particle.species][particle.index];
        }

    /**
     * Returns the vector containing the position of each particle
     */
    const std::vector<unsigned int> & r_vector (unsigned int species) const
        {
            BOOST_ASSERT(species < get_N_species());
            return r[species];
        }

    /**
     * Returns true if the particle is a valid (index, species) combination
     */
    bool particle_is_valid (Particle particle) const
        {
            return (particle.species < get_N_species()
                    && particle.index < get_N_filled(particle.species));
        }

    /** moves a particle to a new position
     *
     * @param particle particle index
     * @param position position index of target site
     */
    void update_position (Particle particle, unsigned int position)
        {
            BOOST_ASSERT(particle_is_valid(particle));
            BOOST_ASSERT(position < get_N_sites());

            const unsigned int old_position = (*this)[particle];

            BOOST_ASSERT(is_occupied(old_position, particle.species));
            // enforce pauli exclusion
            BOOST_ASSERT(!is_occupied(position, particle.species)
                         || r[particle.species][particle.index] == position);

            --positions[particle.species][old_position];
            ++positions[particle.species][position];
            r[particle.species][particle.index] = position;
        }

    /**
     * Returns true if the given position contains a particle of the given
     * species
     */
    bool is_occupied (unsigned int position, unsigned int species) const
        {
            BOOST_ASSERT(position < get_N_sites());
            BOOST_ASSERT(species < get_N_species());
            return bool(positions[species][position] != 0);
        }
    
    /**
     * Returns the index of a given particle of type species located at position;
     * if no such particle is present, returns -1
     */
    int particle_index_at_pos (unsigned int position, unsigned int species) const
        {
            BOOST_ASSERT(position < get_N_sites());
            BOOST_ASSERT(species < get_N_species());

            if (positions[species][position] == 0)
                return -1;
            for (unsigned int i = 0; i < get_N_filled(species); ++i) {
                if (r[species][i] == position)
                    return i;
            }

            // should not be reached
            BOOST_ASSERT(false);
            return -1;
        }

    /**
     * Returns the number of sites on the lattice
     */
    unsigned int get_N_sites (void) const
        {
            BOOST_ASSERT(positions.size() != 0);
            return positions[0].size();
        }

    /**
     * Returns the number of particles of a given species
     */
    unsigned int get_N_filled (unsigned int species) const
        {
            BOOST_ASSERT(species < get_N_species());
            return r[species].size();
        }

    /**
     * Returns the total number of particles of all species
     */
    unsigned int get_N_filled_total (void) const
        {
            return N_filled_total;
        }

    /**
     * Returns the number of species of particles
     */
    unsigned int get_N_species (void) const
        {
            return r.size();
        }

private:
    void _populate_positions (unsigned int N_sites);

    std::vector<std::vector<unsigned int> > r; /**< stores the position of each particle (a vector for each species) */
    std::vector<std::vector<int> > positions; /**< stores the particle(s) on each site (a vector for each species) */

    unsigned int N_filled_total;
};

#endif

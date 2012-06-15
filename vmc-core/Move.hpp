#ifndef _VMC_MOVE_HPP
#define _VMC_MOVE_HPP

#include "lw_vector.hpp"
#include "PositionArguments.hpp"

struct SingleParticleMove
{
    Particle particle;
    unsigned int destination;

    SingleParticleMove (const Particle &particle_, unsigned int destination_)
        : particle(particle_),
          destination(destination_)
        {
        }

    // default constructor is needed so lw_vector can work using boost::array
    SingleParticleMove (void)
        {
        }
};

struct Move : public lw_vector<SingleParticleMove, 2>
{
    bool is_valid_for (const PositionArguments &r) const;
};

#endif

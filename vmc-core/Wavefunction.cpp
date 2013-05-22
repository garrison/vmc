#include <vector>

#include "Wavefunction.hpp"
#include "random-configuration.hpp"
#include "random-move.hpp"

template <typename AmplitudeType>
void Wavefunction<AmplitudeType>::Amplitude::perform_move (const Move &move)
{
    BOOST_ASSERT(!move_in_progress);

    // there's no reason to be making null moves, and if we explicitly disallow
    // them the subclasses can be cleaner.
    BOOST_ASSERT(move.size() != 0);

    BOOST_ASSERT(move.is_valid_for(r));

    // perform the move in the PositionArguments r while remembering the
    // reverse move
    reverse_move = move;
    for (unsigned int i = 0; i < move.size(); ++i) {
        reverse_move[i].destination = r[move[i].particle];
        r.update_position(move[i].particle, move[i].destination);
    }

    BOOST_ASSERT(reverse_move.is_valid_for(r));

    // have the subclass update things such that psi_() returns the new
    // amplitude
    perform_move_(move);

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    move_in_progress = true;
#endif
}

template <typename AmplitudeType>
void Wavefunction<AmplitudeType>::Amplitude::cancel_move (void)
{
    BOOST_ASSERT(move_in_progress);

    // have the subclass cancel the move so psi_() returns the same value as
    // before the attempted move
    cancel_move_();

    // return the PositionArguments r to their original state
    const Move &move = get_reverse_move();
    for (unsigned int i = 0; i < move.size(); ++i)
        r.update_position(move[i].particle, move[i].destination);

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    move_in_progress = false;
#endif
}

template <typename AmplitudeType>
boost::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> Wavefunction<AmplitudeType>::create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const
{
    while (n_attempts--) {
        std::vector<std::vector<unsigned int> > vv;
        for (unsigned int i = 0; i < get_N_species(); ++i)
            vv.push_back(some_random_configuration(get_N_filled(i), *lattice, rng));
        boost::shared_ptr<Wavefunction<AmplitudeType>::Amplitude> wfa(create_wavefunctionamplitude(this_ptr, PositionArguments(vv, lattice->total_sites())));
        if (wfa->is_nonzero())
            return wfa;
    }
    return boost::shared_ptr<Wavefunction<AmplitudeType>::Amplitude>();
}

template <typename AmplitudeType>
Move Wavefunction<AmplitudeType>::Amplitude::propose_move (RandomNumberGenerator &rng) const
{
    // by default we attempt to move a random particle to an empty
    // site. subclasses can feel free to override this behavior, but balance
    // must be maintained!
    Move move;
    const Particle particle(choose_random_particle(r, rng));
    const unsigned int proposed_site_index = plan_particle_move_to_nearby_empty_site(particle, r, *wf->lattice, rng);
    if (proposed_site_index != r[particle])
        move.push_back(SingleParticleMove(particle, proposed_site_index));
    return move;
}

#define VMC_SUPPORTED_TYPE(type) template class Wavefunction<type>
#include "vmc-supported-types.hpp"

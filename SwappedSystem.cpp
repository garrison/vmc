#include "SwappedSystem.hpp"
#include "Lattice.hpp"

// REMEMBER: make sure phibeta's are always updated with copy-on-write

template <typename T>
static int vector_find (const std::vector<T> &vec, const T &val)
{
    for (unsigned int i = 0; i < vec.size(); ++i) {
        if (vec[i] == val)
            return i;
    }
    return -1;
}

SwappedSystem::SwappedSystem (const boost::shared_ptr<const Subsystem> &subsystem_)
    : subsystem(subsystem_),
      phibeta1_dirty(false),
      phibeta2_dirty(false),
      next_step(INITIALIZE)
{
}

void SwappedSystem::initialize (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2)
{
    BOOST_ASSERT(next_step == INITIALIZE);

    const PositionArguments &r1 = phialpha1.get_positions();
    const PositionArguments &r2 = phialpha2.get_positions();

    // FIXME: we need a way to assert that phialpha1 and phialpha2 represent
    // the same wave function, just with different amplitudes.  Then again,
    // we're only calling this function from two places in the code where this
    // can be easily verified ...

    BOOST_ASSERT(r1.get_N_filled() == r2.get_N_filled());
    BOOST_ASSERT(r1.get_N_sites() == r2.get_N_sites());
    BOOST_ASSERT(subsystem->lattice_makes_sense(phialpha1.get_lattice()));
    BOOST_ASSERT(subsystem->lattice_makes_sense(phialpha2.get_lattice()));
    BOOST_ASSERT(&phialpha1.get_lattice() == &phialpha2.get_lattice());

    unsigned int N = r1.get_N_filled();
    for (unsigned int i = 0; i < N; ++i) {
        if (subsystem->particle_is_within(r1[i], phialpha1.get_lattice()))
            copy1_subsystem_indices.push_back(i);
        if (subsystem->particle_is_within(r2[i], phialpha2.get_lattice()))
            copy2_subsystem_indices.push_back(i);
    }

    if (copy1_subsystem_indices.size() == copy2_subsystem_indices.size())
        reinitialize_phibetas(phialpha1, phialpha2);

    next_step = UPDATE;
}

void SwappedSystem::update (int index1, int index2, const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2)
{
    // this function should be called *after* the phialpha's have been updated

    BOOST_ASSERT(next_step == UPDATE);
    next_step = FINISH_UPDATE;

    const PositionArguments &r1 = phialpha1.get_positions();
    const PositionArguments &r2 = phialpha2.get_positions();
    BOOST_ASSERT(r1.size() == r2.size());

    // index1 and index2 are the indices within r1, r2 of the particles that moved
    BOOST_ASSERT(index1 >= -1 && index1 < (int) r1.size());
    BOOST_ASSERT(index2 >= -1 && index2 < (int) r2.size());

    const Lattice &lattice = phialpha1.get_lattice();

    const bool phibetas_were_up_to_date = (get_N_subsystem1() == get_N_subsystem2());

    // these will be will be >= 0 if the particle was in the subsystem before,
    // -1 if the particle was not in the subsystem before, and -2 if the
    // particle isn't even being moved.
    const int pairing_index1 = (index1 == -1) ? -2 : vector_find(copy1_subsystem_indices, (unsigned int) index1);
    const int pairing_index2 = (index2 == -1) ? -2 : vector_find(copy2_subsystem_indices, (unsigned int) index2);

    const bool particle1_now_in_subsystem = (index1 != -1 && subsystem->particle_is_within(r1[index1], lattice));
    const bool particle2_now_in_subsystem = (index2 != -1 && subsystem->particle_is_within(r2[index2], lattice));

    const int delta1 = (particle1_now_in_subsystem ? 1 : 0) + (pairing_index1 >= 0 ? -1 : 0);
    const int delta2 = (particle2_now_in_subsystem ? 1 : 0) + (pairing_index2 >= 0 ? -1 : 0);

    if (index1 == -1)
        BOOST_ASSERT(delta1 == 0);
    if (index2 == -1)
        BOOST_ASSERT(delta2 == 0);

    if (phibetas_were_up_to_date && delta1 == -1 && delta2 == -1) {
        // if the phibetas are up to date and a particle leaves each subsystem
        // simultaneously, we need to use some special logic in case we have to
        // re-pair the remaining particles in the subsystem.  (re-pair in the
        // sense of what gets swapped with what)

        // these are equivalent to the "if" statement, but are a useful sanity
        // check nonetheless (for now)
        BOOST_ASSERT(pairing_index1 >= 0 && pairing_index2 >= 0);
        BOOST_ASSERT(!particle1_now_in_subsystem);
        BOOST_ASSERT(!particle2_now_in_subsystem);

        // copy-on-write for the phibeta's
        if (!phibeta1.unique())
            phibeta1 = phibeta1->clone();
        if (!phibeta2.unique())
            phibeta2 = phibeta2->clone();

        // update the phibeta's for the particles that left the system
        BOOST_ASSERT(!phibeta1_dirty && !phibeta2_dirty);
        phibeta1->move_particle(index1, r1[index1]);
        phibeta2->move_particle(index2, r2[index2]);
        phibeta1_dirty = true;
        phibeta2_dirty = true;

        if (pairing_index1 != pairing_index2) { // we must re-pair
            // update the phibeta's for the particles that are about to be paired
            BOOST_ASSERT(phibeta1_dirty && phibeta2_dirty);
            phibeta1->finish_particle_moved_update();
            phibeta2->finish_particle_moved_update();
            phibeta1->move_particle(copy1_subsystem_indices[pairing_index2], r2[copy2_subsystem_indices[pairing_index1]]);
            phibeta2->move_particle(copy2_subsystem_indices[pairing_index1], r1[copy1_subsystem_indices[pairing_index2]]);

            // update the subsystem indices so they become paired at the min_pairing_index
            if (pairing_index1 < pairing_index2)
                copy1_subsystem_indices[pairing_index1] = copy1_subsystem_indices[pairing_index2];
            else
                copy2_subsystem_indices[pairing_index2] = copy2_subsystem_indices[pairing_index1];
        }

        // remove the empty pair in the subsystem indices (yes, these steps
        // make sense whether we had to re-pair or not)
        const unsigned int max_pairing_index = std::max(pairing_index1, pairing_index2);
        copy1_subsystem_indices[max_pairing_index] = copy1_subsystem_indices[copy1_subsystem_indices.size() - 1];
        copy2_subsystem_indices[max_pairing_index] = copy2_subsystem_indices[copy2_subsystem_indices.size() - 1];
        copy1_subsystem_indices.pop_back();
        copy2_subsystem_indices.pop_back();
    } else {
        // update the subsystem indices
        if (delta1 == 1) {
            copy1_subsystem_indices.push_back((unsigned int) index1);
        } else if (delta1 == -1) {
            copy1_subsystem_indices[pairing_index1] = copy1_subsystem_indices[copy1_subsystem_indices.size() - 1];
            copy1_subsystem_indices.pop_back();
        }
        if (delta2 == 1) {
            copy2_subsystem_indices.push_back((unsigned int) index2);
        } else if (delta2 == -1) {
            copy2_subsystem_indices[pairing_index2] = copy2_subsystem_indices[copy2_subsystem_indices.size() - 1];
            copy2_subsystem_indices.pop_back();
        }

        // if there is a different number of particles in each subsystem, we
        // can't update the phibeta's, so we forget them for now and return.
        if (get_N_subsystem1() != get_N_subsystem2()) {
            phibeta1.reset();
            phibeta2.reset();
            return;
        }

        if (!phibetas_were_up_to_date) {
            reinitialize_phibetas(phialpha1, phialpha2);
        } else {
            // either both particles moved within their respective subsystems
            // (if they moved at all), or both entered the subsystem and paired
            // with each other immediately.
            BOOST_ASSERT(delta1 == delta2);
            BOOST_ASSERT(delta1 == 0 || delta1 == 1);

            // update the phibeta's, performing copy-on-write
            if (index1 != -1) {
                boost::shared_ptr<WavefunctionAmplitude> &phibeta = particle1_now_in_subsystem ? phibeta2 : phibeta1;
                bool &phibeta_dirty = particle1_now_in_subsystem ? phibeta2_dirty : phibeta1_dirty;
                const unsigned int phibeta_index = particle1_now_in_subsystem ? copy2_subsystem_indices[(delta1 == 0) ? (unsigned int) pairing_index1 : copy1_subsystem_indices.size() - 1] : (unsigned int) index1;
                if (!phibeta.unique())
                    phibeta = phibeta->clone();
                BOOST_ASSERT(!phibeta_dirty); // will always be clean here, but not necessarily below
                phibeta->move_particle(phibeta_index, r1[index1]);
                phibeta_dirty = true;
            }

            if (index2 != -1) {
                boost::shared_ptr<WavefunctionAmplitude> &phibeta = particle2_now_in_subsystem ? phibeta1 : phibeta2;
                bool &phibeta_dirty = particle2_now_in_subsystem ? phibeta1_dirty : phibeta2_dirty;
                const unsigned int phibeta_index = particle2_now_in_subsystem ? copy1_subsystem_indices[(delta2 == 0) ? (unsigned int) pairing_index2 : copy2_subsystem_indices.size() - 1] : (unsigned int) index2;
                if (!phibeta.unique())
                    phibeta = phibeta->clone();
                if (phibeta_dirty)
                    phibeta->finish_particle_moved_update();
                phibeta->move_particle(phibeta_index, r2[index2]);
                phibeta_dirty = true;
            }
        }
    }
}

void SwappedSystem::finish_update (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2)
{
    BOOST_ASSERT(next_step == FINISH_UPDATE);
    next_step = UPDATE;

    if (phibeta1_dirty) {
        BOOST_ASSERT(get_N_subsystem1() == get_N_subsystem2());
        phibeta1->finish_particle_moved_update();
    }
    phibeta1_dirty = false;

    if (phibeta2_dirty) {
        BOOST_ASSERT(get_N_subsystem1() == get_N_subsystem2());
        phibeta2->finish_particle_moved_update();
    }
    phibeta2_dirty = false;

#ifdef CAREFUL
    if (get_N_subsystem1() == get_N_subsystem2())
        verify_phibetas(phialpha1, phialpha2);
#else
    (void) phialpha1;
    (void) phialpha2;
#endif
}

void SwappedSystem::reinitialize_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2)
{
    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());

#ifdef DEBUG
    std::cerr << "swapping " << copy1_subsystem_indices.size() << " particles" << std::endl;
#endif

    PositionArguments swapped_r1(phialpha1.get_positions()), swapped_r2(phialpha2.get_positions());
    swap_positions(swapped_r1, swapped_r2);

    phibeta1 = phialpha1.clone();
    phibeta1->reset(swapped_r1);
    phibeta1_dirty = false;

    phibeta2 = phialpha2.clone();
    phibeta2->reset(swapped_r2);
    phibeta2_dirty = false;

#ifdef CAREFUL
    verify_phibetas(phialpha1, phialpha2);
#endif
}

void SwappedSystem::verify_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2) const
{
#ifdef BOOST_DISABLE_ASSERTS
    (void) phialpha1;
    (void) phialpha2;
#else
    const PositionArguments &r1 = phialpha1.get_positions();
    const PositionArguments &r2 = phialpha2.get_positions();

    BOOST_ASSERT(r1.get_N_filled() == r2.get_N_filled());

    const unsigned int N = r1.get_N_filled();
    const Lattice &lattice = phialpha1.get_lattice();

    // verify that the system index arrays have everything they need (and no duplicates!)
    unsigned int c1 = 0, c2 = 0;
    for (unsigned int i = 0; i < N; ++i) {
        const bool b1 = vector_find(copy1_subsystem_indices, i) != -1;
        const bool b2 = vector_find(copy2_subsystem_indices, i) != -1;
        if (b1) ++c1;
        if (b2) ++c2;
        BOOST_ASSERT(b1 == subsystem->particle_is_within(r1[i], lattice));
        BOOST_ASSERT(b2 == subsystem->particle_is_within(r2[i], lattice));
    }
    BOOST_ASSERT(c1 == copy1_subsystem_indices.size());
    BOOST_ASSERT(c2 == copy2_subsystem_indices.size());

    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());
    BOOST_ASSERT(phibeta1 != 0);
    BOOST_ASSERT(phibeta2 != 0);

    // verify that the positions in the phibeta's are correct
    PositionArguments swapped_r1(phialpha1.get_positions()), swapped_r2(phialpha2.get_positions());
    swap_positions(swapped_r1, swapped_r2);

    for (unsigned int i = 0; i < N; ++i) {
        BOOST_ASSERT(swapped_r1[i] == phibeta1->get_positions()[i]);
        BOOST_ASSERT(swapped_r2[i] == phibeta2->get_positions()[i]);
    }
#endif
}

// this is a private utility method, called from both reinitialize_phibetas()
// and verify_phibetas()
void SwappedSystem::swap_positions (PositionArguments &r1, PositionArguments &r2) const
{
    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());

    const unsigned int particles_to_swap = copy1_subsystem_indices.size();
    for (unsigned int i = 0; i < particles_to_swap; ++i) {
        unsigned int i1 = copy1_subsystem_indices[i], i2 = copy2_subsystem_indices[i];
        // perform the swap
        unsigned int tmp_pos = r1[i1];
        r1.update_position(i1, r2[i2]);
        r2.update_position(i2, tmp_pos);
    }
}

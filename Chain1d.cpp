#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "Chain1d.hpp"
#include "random-combination.hpp"

Chain1dArguments::Chain1dArguments (const Chain1d &wf_, rng_class &rng)
{
    random_combination(r, wf_.get_N_filled(), wf_.get_N_sites(), rng);
}

Chain1dArguments::Chain1dArguments (const Chain1d &wf_, const std::vector<Chain1d::position_t> &r_)
    : r(r_)
{
    BOOST_ASSERT((int)r_.size() == wf_.get_N_filled());
    // fixme: assert each position is valid
}

static int move_random_particle_randomly(Chain1dArguments &r, const Chain1d &wf, rng_class &rng, bool adjacent_only=true)
{
    boost::uniform_smallint<> integer_distribution(0, wf.get_N_filled() - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > particle_gen(rng, integer_distribution);
    boost::uniform_smallint<> movement_distribution(0, adjacent_only ? 1 : (wf.get_N_sites() - 1));
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > movement_gen(rng, movement_distribution);

    // choose particle
    int chosen_particle = particle_gen();
#ifdef DEBUG
    std::cerr << "Moving particle " << chosen_particle << " from " << r[chosen_particle];
#endif

    // move it randomly
    Chain1d::position_t rcp;
    if (adjacent_only) {
	rcp = r[chosen_particle];
	rcp += movement_gen() * 2 - 1;
	// enforce PBC
	if (rcp >= wf.get_N_sites())
	    rcp -= wf.get_N_sites();
	else if (rcp < 0)
	    rcp += wf.get_N_sites();
    } else {
	// fixme: choose any random *empty* site
	rcp = movement_gen();
    }
    r.update_position(chosen_particle, rcp);
#ifdef DEBUG
    std::cerr << " to " << rcp << std::endl;
#endif

    return chosen_particle;
}

Chain1dWalk::Chain1dWalk (const Chain1d &wf_, const Chain1dArguments &arguments_)
    : wf(new Chain1d(wf_)),
      r(arguments_),
      transition_in_progress(false)
{
    //BOOST_ASSERT(wf_ == r.wf);
    int N = wf->get_N_filled();
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j)
	    mat(i, j) = wf->phi(j, r[i]);
    }
    cmat = mat;
}

probability_t Chain1dWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_in_progress);

    int N = wf->get_N_filled();

    int chosen_particle = move_random_particle_randomly(r, *wf, rng);

    // calculate each phi at new position and update Ceperly matrix
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, 1> phivec(N);
    for (int i = 0; i < N; ++i)
	phivec(i) = wf->phi(i, r[chosen_particle]);
    cmat.update_row(chosen_particle, phivec);

    // take ratio of determinants and return a probability
    Chain1d::amplitude_t amplitude_ratio = cmat.calculate_determinant_ratio();
    probability_t rv = norm(amplitude_ratio);
#ifdef DEBUG
    std::cerr << "ratio " << rv << std::endl;
#endif
    transition_in_progress = true;
    return rv;
}

void Chain1dWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

#ifdef DEBUG
    for (unsigned int i = 0; i < r.size(); ++i)
	std::cerr << r[i] << ' ';
    std::cerr << std::endl;
#endif

    cmat.finish_row_update();

    transition_in_progress = false;
}

Chain1dRenyiWalk::Chain1dRenyiWalk (const Chain1d &wf_, const Subsystem<Chain1d> *subsystem_, RenyiWalkType walk_type_, rng_class &rng)
    : wf(new Chain1d(wf_)),
      r1(wf_, rng),
      r2(r1),
      subsystem(subsystem_),
      N_subsystem1(fermion_partition<Chain1d>(r1, *subsystem).first),
      N_subsystem2(fermion_partition<Chain1d>(r2, *subsystem).first),
      walk_type(walk_type_),
      transition_in_progress(false)
{
    BOOST_ASSERT(N_subsystem1 == N_subsystem2);

    int N = wf->get_N_filled();
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat1(N, N), mat2(N, N);
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    mat1(i, j) = wf->phi(j, r1[i]);
	    mat2(i, j) = wf->phi(j, r2[i]);
	}
    }
    phialpha1 = mat1;
    phialpha2 = mat2;

    for (int i = 0; i < N_subsystem1; ++i)
	mat1.row(i).swap(mat2.row(i));
    phibeta1 = mat1;
    phibeta2 = mat2;
}

static void consider_crossing (const Subsystem<Chain1d> *subsystem, Chain1dArguments &r1,
			       int &chosen_particle1, int &N_subsystem1,
			       CeperlyMatrix<Chain1d::amplitude_t> &phialpha1,
			       CeperlyMatrix<Chain1d::amplitude_t> &phibeta1,
			       CeperlyMatrix<Chain1d::amplitude_t> &phibeta2)
{
    if (chosen_particle1 < N_subsystem1) {
	// chosen particle was in subsystem at beginning of step
	if (!subsystem->particle_is_within(r1[chosen_particle1])) {
	    // particle moved out of subsystem
#ifdef DEBUG
	    std::cerr << "Particle leaving subsystem" << std::endl;
#endif
	    --N_subsystem1;
	    if (chosen_particle1 != N_subsystem1) {
		// swap two particles that were inside the subsystem so that
		// the chosen_particle is on the border
		r1.swap_positions(chosen_particle1, N_subsystem1);
		phialpha1.swap_rows(chosen_particle1, N_subsystem1);
		phibeta2.swap_rows(chosen_particle1, N_subsystem1);
		chosen_particle1 = N_subsystem1;
	    }
	    //return true;
	}
    } else {
	// chosen particle was not in subsystem at beginning of step
	if (subsystem->particle_is_within(r1[chosen_particle1])) {
	    // particle entered subsystem
#ifdef DEBUG
	    std::cerr << "Particle entering subsystem" << std::endl;
#endif
	    if (chosen_particle1 != N_subsystem1) {
		// swap two particles that were outside the subsystem so that
		// the chosen_particle is on the border
		r1.swap_positions(chosen_particle1, N_subsystem1);
		phialpha1.swap_rows(chosen_particle1, N_subsystem1);
		phibeta1.swap_rows(chosen_particle1, N_subsystem1);
		chosen_particle1 = N_subsystem1;
	    }
	    ++N_subsystem1;
	    //return true;
	}
    }
    //return false;
}

probability_t Chain1dRenyiWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_in_progress);
    transition_in_progress = true;

    int N = wf->get_N_filled();

    boost::uniform_smallint<> adjacent_distribution(0, 30);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > adjacent_gen(rng, adjacent_distribution);
    bool adjacent_only = (adjacent_gen() > 0);

    int chosen_particle1 = move_random_particle_randomly(r1, *wf, rng, adjacent_only);
    int chosen_particle2 = move_random_particle_randomly(r2, *wf, rng, adjacent_only);

    // if a particle crossed the subsystem boundary, account for that
    consider_crossing(subsystem, r1, chosen_particle1, N_subsystem1,
		      phialpha1, phibeta1, phibeta2);
    consider_crossing(subsystem, r2, chosen_particle2, N_subsystem2,
		      phialpha2, phibeta2, phibeta1);

    // automatic reject if the subsystems now have different particle counts
    if (N_subsystem1 != N_subsystem2)
	return 0;

    // calculate each phi at new position and update phialpha Ceperly matrices
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, 1> phivec1(N), phivec2(N);

    for (int i = 0; i < N; ++i) {
	phivec1(i) = wf->phi(i, r1[chosen_particle1]);
	phivec2(i) = wf->phi(i, r2[chosen_particle2]);
    }
    phialpha1.update_row(chosen_particle1, phivec1);
    phialpha2.update_row(chosen_particle2, phivec2);

    // update phibeta matrices
    //
    // fixme: this is repetitive and error-prone, so we could probably abstract
    // this away a little bit ...
    Chain1d::amplitude_t phibeta1_ratio = 1, phibeta2_ratio = 1;
    if (chosen_particle1 < N_subsystem1) {
	phibeta2.update_row(chosen_particle1, phivec1);
	phibeta2_ratio *= phibeta2.calculate_determinant_ratio();
	phibeta2.finish_row_update();
    } else {
	phibeta1.update_row(chosen_particle1, phivec1);
	phibeta1_ratio *= phibeta1.calculate_determinant_ratio();
	phibeta1.finish_row_update();
    }
    if (chosen_particle2 < N_subsystem2) {
	phibeta1.update_row(chosen_particle2, phivec2);
	phibeta1_ratio *= phibeta1.calculate_determinant_ratio();
	phibeta1.finish_row_update();
    } else {
	phibeta2.update_row(chosen_particle2, phivec2);
	phibeta2_ratio *= phibeta2.calculate_determinant_ratio();
	phibeta2.finish_row_update();
    }

    // calculate ratio of determinants and return a probability
    if (walk_type == SWAPA_MOD) {
	return (norm(phialpha1.calculate_determinant_ratio()) *
		norm(phialpha2.calculate_determinant_ratio()));
    } else { // walk_type == SWAPA_SIGN
	return abs(phialpha1.calculate_determinant_ratio() *
		   phialpha2.calculate_determinant_ratio() *
		   phibeta1_ratio *
		   phibeta2_ratio);
    }
}

void Chain1dRenyiWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

#ifdef DEBUG
    for (unsigned int i = 0; i < r1.size(); ++i)
	std::cerr << r1[i] << ' ';
    std::cerr << std::endl;
    for (unsigned int i = 0; i < r2.size(); ++i)
	std::cerr << r2[i] << ' ';
    std::cerr << std::endl << std::endl;
#endif

    phialpha1.finish_row_update();
    phialpha2.finish_row_update();

    transition_in_progress = false;
}

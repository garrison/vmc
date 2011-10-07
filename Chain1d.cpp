#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "Chain1d.hpp"

Chain1dWalk::Chain1dWalk (const Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> &mat, const Chain1d &wf_, const std::vector<int> &r_)
    : wf(new Chain1d(wf_)),
      r(r_),
      cmat(mat),
      transition_in_progress(false)
{
    BOOST_ASSERT((int)r.size() == wf->get_N_filled());
}

std::auto_ptr<Chain1dWalk> Chain1dWalk::random_initial_state (const Chain1d &wf_)
{
    int N = wf_.get_N_filled();

    //boost::mt19937 rng(0); // fixme: seed
    //boost::uniform_01<boost::mt19937> zeroone(rng);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    std::vector<int> r(N); // fixme
    for (int i = 0; i < N; ++i) {
	r[i] = i; //zeroone();
	for (int j = 0; j < N; ++j)
	    mat(i, j) = wf_.phi(j, r[i]);
    }

    return std::auto_ptr<Chain1dWalk>(new Chain1dWalk(mat, wf_, r));
}

probability_t Chain1dWalk::compute_probability_ratio_of_random_transition (void)
{
    BOOST_ASSERT(!transition_in_progress);

    int N = wf->get_N_filled();
    static boost::mt19937 generator(4); // fixme: seed
    static boost::uniform_smallint<> integer_distribution(0, N - 1);
    static boost::variate_generator<boost::mt19937&, boost::uniform_smallint<> > particle_gen(generator, integer_distribution);
    static boost::uniform_smallint<> binary_distribution(0, 1);
    static boost::variate_generator<boost::mt19937&, boost::uniform_smallint<> > movement_gen(generator, binary_distribution);
    //static boost::uniform_real<> movement_distribution(-0.1, 0.1);
    //static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > distance_gen(generator, movement_distribution);

    int chosen_particle = particle_gen();
    int *rcp = &r[chosen_particle];
    std::cerr << "Moving particle " << chosen_particle << " from " << *rcp;

    // move it randomly
    *rcp += movement_gen() * 2 - 1;
    // enforce PBC
    if (*rcp >= wf->get_N_sites())
	*rcp -= wf->get_N_sites();
    else if (*rcp < 0)
	*rcp += wf->get_N_sites();
    std::cerr << " to " << *rcp << std::endl;

    // calculate each phi at new position
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, 1> phivec(N);
    for (int i = 0; i < N; ++i)
	phivec(i) = wf->phi(i, *rcp);

    cmat.update_row(chosen_particle, phivec);

    // take ratio of determinants and return a probability
    amplitude_t amplitude_ratio = cmat.calculate_determinant_ratio();
    probability_t rv = norm(amplitude_ratio);
    std::cerr << "ratio " << rv << std::endl;
    transition_in_progress = true;
    return rv;
}

void Chain1dWalk::finalize_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    for (int i = 0; i < r.size(); ++i)
	std::cerr << r[i] << ' ';
    std::cerr << std::endl;

    cmat.finish_row_update();

    transition_in_progress = false;
}

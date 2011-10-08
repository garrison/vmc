#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "Chain1d.hpp"

Chain1dWalk::Chain1dWalk (const Chain1d &wf_)
    : wf(new Chain1d(wf_)),
      r(wf->get_N_filled()),
      transition_in_progress(false)
{
    for (int i = 0; i < wf->get_N_filled(); ++i)
	r[i] = i; // fixme!

    initialize_cmat();
}

Chain1dWalk::Chain1dWalk (const Chain1d &wf_, const std::vector<Chain1d::position_t> &r_)
    : wf(new Chain1d(wf_)),
      r(r_),
      transition_in_progress(false)
{
    BOOST_ASSERT((int)r.size() == wf->get_N_filled());
    initialize_cmat();
}

void Chain1dWalk::initialize_cmat (void)
{
    int N = wf->get_N_filled();
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j)
	    mat(i, j) = wf->phi(j, r[i]);
    }
    cmat = mat;
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
    std::cerr << "Moving particle " << chosen_particle << " from " << r[chosen_particle];

    // move it randomly
    r[chosen_particle] += movement_gen() * 2 - 1;
    // enforce PBC
    if (r[chosen_particle] >= wf->get_N_sites())
	r[chosen_particle] -= wf->get_N_sites();
    else if (r[chosen_particle] < 0)
	r[chosen_particle] += wf->get_N_sites();
    std::cerr << " to " << r[chosen_particle] << std::endl;

    // calculate each phi at new position
    Eigen::Matrix<Chain1d::amplitude_t, Eigen::Dynamic, 1> phivec(N);
    for (int i = 0; i < N; ++i)
	phivec(i) = wf->phi(i, r[chosen_particle]);

    cmat.update_row(chosen_particle, phivec);

    // take ratio of determinants and return a probability
    Chain1d::amplitude_t amplitude_ratio = cmat.calculate_determinant_ratio();
    probability_t rv = norm(amplitude_ratio);
    std::cerr << "ratio " << rv << std::endl;
    transition_in_progress = true;
    return rv;
}

void Chain1dWalk::finalize_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    for (unsigned int i = 0; i < r.size(); ++i)
	std::cerr << r[i] << ' ';
    std::cerr << std::endl;

    cmat.finish_row_update();

    transition_in_progress = false;
}

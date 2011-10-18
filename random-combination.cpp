#include <algorithm>
#include <set>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "random-combination.hpp"
#include "vmc-typedefs.hpp"

// http://stackoverflow.com/questions/2394246/algorithm-to-select-a-single-random-combination-of-values
void random_combination (std::vector<int> &v, int r, int n, rng_class &rng)
{
    // per Jon Bentley's article in CACM, September 1987, Volume 30, Number 9
    BOOST_ASSERT(n > 0);
    BOOST_ASSERT(r > 0);
    BOOST_ASSERT(r <= n);

    std::set<int> vs;
    v.resize(0);
    v.reserve(r);

    for (int k = n - r; k < n; ++k) {
	boost::uniform_smallint<> dist(0, k - 1);
	// fixme: make sure we can use the same rng again and again here without updating it
	boost::variate_generator<rng_class&, boost::uniform_smallint<> > gen(rng, dist);
	int x = gen();
	int a = (vs.find(x) != vs.end()) ? k : x;
	v.push_back(a);
	vs.insert(a);
    }

    std::sort_heap(v.begin(), v.end());
}


#include <boost/assert.hpp>

#include "JordanWignerJastrowFactor.hpp"

template <typename AmplitudeType>
Big<AmplitudeType> JordanWignerJastrowFactor<AmplitudeType>::compute_jastrow (const PositionArguments &r) const
{
    // FIXME: we should relax this
    BOOST_ASSERT(r.get_N_species() == 1);

    unsigned int c = 0;
    for (unsigned int i = 0; i < r.get_N_filled(0); ++i) {
        const Particle pi(i, 0);
        for (unsigned int j = i + 1; j < r.get_N_filled(0); ++j) {
            const Particle pj(j, 0);
            BOOST_ASSERT(r[pi] != r[pj]);
            if (r[pi] < r[pj])
                ++c;
        }
    }
    const int sgn = (c & 1) ? -1 : 1;
    return Big<AmplitudeType>(sgn);
}

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) template class JordanWignerJastrowFactor<type>
#include "vmc-supported-types.hpp"

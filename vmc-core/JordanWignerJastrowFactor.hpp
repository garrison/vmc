#ifndef _VMC_JORDAN_WIGNER_JASTROW_FACTOR
#define _VMC_JORDAN_WIGNER_JASTROW_FACTOR

#include "JastrowFactor.hpp"

/**
 * Performs Jordan-Wigner in first quantization
 *
 * Currently, should only be used when get_N_species == 1.
 */
template <typename AmplitudeType>
class JordanWignerJastrowFactor : public JastrowFactor<AmplitudeType>
{
    virtual Big<AmplitudeType> compute_jastrow (const PositionArguments &r) const override;
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class JordanWignerJastrowFactor<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif

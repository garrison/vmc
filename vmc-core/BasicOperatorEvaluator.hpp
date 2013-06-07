#ifndef _VMC_BASIC_OPERATOR_EVALUATOR_HPP
#define _VMC_BASIC_OPERATOR_EVALUATOR_HPP

#include "BoundaryCondition.hpp"
#include "BasicOperator.hpp"
#include "Wavefunction.hpp"

template <typename _AmplitudeType>
class BasicOperatorEvaluator
{
public:
    typedef _AmplitudeType AmplitudeType;
    typedef AmplitudeType PhaseType;

    const BasicOperator m_operator;
    const BoundaryConditions<PhaseType> bcs;

    BasicOperatorEvaluator (const BasicOperator &operator_, const BoundaryConditions<PhaseType> &bcs_);

    AmplitudeType evaluate (const typename Wavefunction<AmplitudeType>::Amplitude &wfa) const;

private:
    bool is_sum_over_sites (void) const
        {
            return bcs.size() != 0;
        }

    const unsigned int min_required_species;
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class BasicOperatorEvaluator<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif

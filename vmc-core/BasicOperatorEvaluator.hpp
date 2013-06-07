#ifndef _VMC_BASIC_OPERATOR_EVALUATOR_HPP
#define _VMC_BASIC_OPERATOR_EVALUATOR_HPP

#include "BoundaryCondition.hpp"
#include "BasicOperator.hpp"
#include "Wavefunction.hpp"

class BasicOperatorEvaluator
{
public:
    const BasicOperator m_operator;

    BasicOperatorEvaluator (const BasicOperator &operator_);

    /**
     * Evaluate the operator
     *
     * If boundary conditions are given, this will loop over all sites, summing
     * the expectation value of the operator as translated across the lattice
     * (wrapping around any boundary conditions that are not open).  If
     * `boundary_conditions` is the empty array, no sum is performed.
     */
    template <typename AmplitudeType>
    AmplitudeType evaluate (const typename Wavefunction<AmplitudeType>::Amplitude &wfa, const BoundaryConditions<AmplitudeType> &boundary_conditions) const;

private:
    const unsigned int min_required_species;
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template type BasicOperatorEvaluator::evaluate(const typename Wavefunction<type>::Amplitude &wfa, const BoundaryConditions<type> &boundary_conditions) const;
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif

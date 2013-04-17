#include <string>
#include <boost/lexical_cast.hpp>

#include "CeperleyMatrix.hpp"

static inline std::string construct_what_string (const real_t &inverse_error, unsigned int n_smw_updates, const real_t &smallest_detrat)
{
    return (std::string("Unrecoverable large inverse error: ")
            + boost::lexical_cast<std::string>(inverse_error)
            + std::string(" after ")
            + boost::lexical_cast<std::string>(n_smw_updates) + std::string(" updates.")
            + std::string("  smallest detrat: ")
            + boost::lexical_cast<std::string>(smallest_detrat));
}

unrecoverable_matrix_inverse_error::unrecoverable_matrix_inverse_error (real_t inverse_error_, unsigned int n_smw_updates_, real_t smallest_detrat_)
    : std::runtime_error(construct_what_string(inverse_error_, n_smw_updates_, smallest_detrat_)),
      inverse_error(inverse_error_),
      n_smw_updates(n_smw_updates_),
      smallest_detrat(smallest_detrat_)
{
}

template<>
const double CeperleyMatrix<std::complex<double> >::ceperley_determinant_base_upper_cutoff = 1e25;
template<>
const long double CeperleyMatrix<std::complex<long double> >::ceperley_determinant_base_upper_cutoff = 1e25;

template<>
const double CeperleyMatrix<std::complex<double> >::ceperley_determinant_base_lower_cutoff = 1e-25;
template<>
const long double CeperleyMatrix<std::complex<long double> >::ceperley_determinant_base_lower_cutoff = 1e-25;

// if this is set too low, we may not be able to reliably recognize singular
// matrices.  In particular, an abs(detrat) as high as 1.27765e-05 has been
// known to cause problems on the DMetal 48x2 "presentation point"
template<>
const double CeperleyMatrix<std::complex<double> >::ceperley_detrat_lower_cutoff = 1e-4;
template<>
const long double CeperleyMatrix<std::complex<long double> >::ceperley_detrat_lower_cutoff = 1e-4;

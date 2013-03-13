#include "CeperleyMatrix.hpp"

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

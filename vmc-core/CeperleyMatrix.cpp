#include "CeperleyMatrix.hpp"

template<>
const double CeperleyMatrix<std::complex<double> >::ceperley_determinant_upper_cutoff = 1e50;
template<>
const long double CeperleyMatrix<std::complex<long double> >::ceperley_determinant_upper_cutoff = 1e50;

// if this is set too low, we may not be able to reliably recognize singular matrices
template<>
const double CeperleyMatrix<std::complex<double> >::ceperley_determinant_lower_cutoff = 1e-4;
template<>
const long double CeperleyMatrix<std::complex<long double> >::ceperley_determinant_lower_cutoff = 1e-4;

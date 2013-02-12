#include "CeperleyMatrix.hpp"

template<>
const std::complex<double> CeperleyMatrix<std::complex<double> >::ceperley_determinant_cutoff = 1e-6;

template<>
const std::complex<long double> CeperleyMatrix<std::complex<long double> >::ceperley_determinant_cutoff = 1e-6;

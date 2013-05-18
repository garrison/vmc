#include <limits>

#include <boost/version.hpp>

// we always want to include vmc-typedefs.hpp before including Eigen
#include "vmc-typedefs.hpp"
#include <Eigen/Core>

#include "RunInformation.hpp"

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#ifdef __VERSION__
const char * RunInformation::compiler = __VERSION__;
#else
const char * RunInformation::compiler = 0;
#endif

int RunInformation::Precision::digits = std::numeric_limits<real_t>::digits;
int RunInformation::Precision::min_exponent = std::numeric_limits<real_t>::min_exponent;
int RunInformation::Precision::max_exponent = std::numeric_limits<real_t>::max_exponent;

const char * RunInformation::boost_version = BOOST_LIB_VERSION;
const char * RunInformation::eigen_version = STRINGIFY(EIGEN_WORLD_VERSION) "." STRINGIFY(EIGEN_MAJOR_VERSION) "." STRINGIFY(EIGEN_MINOR_VERSION);

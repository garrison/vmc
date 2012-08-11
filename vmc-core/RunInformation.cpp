#include <limits>

#include <boost/version.hpp>
#include <Eigen/Core>

#include "RunInformation.hpp"
#include "vmc-typedefs.hpp"

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

Json::Value RunInformation::json_info (void)
{
    Json::Value rv(Json::objectValue);

#ifdef __VERSION__
    rv["compiler"] = __VERSION__;
#else
    rv["compiler"] = Json::Value(Json::nullValue);
#endif
#if 0
    rv["compiler_commandline"] = STRINGIFY(VMC_COMPILER_COMMANDLINE);
#endif

    rv["precision"]["digits"] = std::numeric_limits<real_t>::digits;
    rv["precision"]["min_exponent"] = std::numeric_limits<real_t>::min_exponent;
    rv["precision"]["max_exponent"] = std::numeric_limits<real_t>::max_exponent;

    rv["version"]["boost"] = BOOST_LIB_VERSION;
    rv["version"]["eigen"] = STRINGIFY(EIGEN_WORLD_VERSION) "." STRINGIFY(EIGEN_MAJOR_VERSION) "." STRINGIFY(EIGEN_MINOR_VERSION);

    return rv;
}

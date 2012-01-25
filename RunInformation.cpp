#include <sys/utsname.h>
#include <limits>

#include <boost/version.hpp>
#include <Eigen/Core>

#include "RunInformation.hpp"
#include "vmc-typedefs.hpp"

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

static Json::Value get_json_uname()
{
    Json::Value rv(Json::objectValue);
    struct utsname buf;
    if (uname(&buf) == 0) {
        rv["sysname"] = buf.sysname;
        rv["nodename"] = buf.nodename;
        rv["release"] = buf.release;
        rv["version"] = buf.version;
        rv["machine"] = buf.machine;
    }
    return rv;
}

static const Json::Value platform = get_json_uname();

RunInformation::RunInformation (void)
{
}

Json::Value RunInformation::json_info (void) const
{
    Json::Value rv(Json::objectValue);
    rv["version"]["boost"] = BOOST_LIB_VERSION;
    rv["version"]["eigen"] = STRINGIFY(EIGEN_WORLD_VERSION) "." STRINGIFY(EIGEN_MAJOR_VERSION) "." STRINGIFY(EIGEN_MINOR_VERSION);
    rv["version"]["vmc"] = "0";
    rv["platform"] = platform;
#ifdef __VERSION__
    rv["compiler"] = __VERSION__;
#else
    rv["compiler"] = Json::Value(Json::nullValue);
#endif
    // fixme: rv["compiler_defines"];
    rv["precision"]["digits"] = std::numeric_limits<real_t>::digits;
    rv["precision"]["min_exponent"] = std::numeric_limits<real_t>::min_exponent;
    rv["precision"]["max_exponent"] = std::numeric_limits<real_t>::max_exponent;
    // fixme: start time, finish time, rusage

    return rv;
}

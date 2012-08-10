#include <sys/utsname.h>
#include <sys/resource.h>
#include <limits>
#include <ctime>
#include <string>

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

static double timeval_to_double (const struct timeval &tv)
{
    return double(tv.tv_sec) + double(.000001 * tv.tv_usec);
}

static std::string time_to_string (const std::time_t &t)
{
    std::string rv(std::ctime(&t));
    rv.erase(rv.size() - 1); // remove the trailing newline
    return rv;
}

RunInformation::RunInformation (void)
{
    std::time(&start_time);
}

Json::Value RunInformation::json_info (void) const
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

    rv["platform"] = platform;

    rv["precision"]["digits"] = std::numeric_limits<real_t>::digits;
    rv["precision"]["min_exponent"] = std::numeric_limits<real_t>::min_exponent;
    rv["precision"]["max_exponent"] = std::numeric_limits<real_t>::max_exponent;

    rv["version"]["boost"] = BOOST_LIB_VERSION;
    rv["version"]["eigen"] = STRINGIFY(EIGEN_WORLD_VERSION) "." STRINGIFY(EIGEN_MAJOR_VERSION) "." STRINGIFY(EIGEN_MINOR_VERSION);

    std::time_t finish_time;
    std::time(&finish_time);
    rv["datetime"]["start"] = time_to_string(start_time);
    rv["datetime"]["finish"] = time_to_string(finish_time);
    rv["walltime"] = Json::Int(finish_time - start_time);

#if 0
    struct rusage res_usage;
    if (getrusage(RUSAGE_SELF, &res_usage) == 0) {
        rv["rusage"]["utime"] = timeval_to_double(res_usage.ru_utime);
        rv["rusage"]["stime"] = timeval_to_double(res_usage.ru_stime);
        rv["rusage"]["maxrss"] = Json::Int(res_usage.ru_maxrss);
        rv["rusage"]["nvcsw"] = Json::Int(res_usage.ru_nvcsw);
        rv["rusage"]["nivcsw"] = Json::Int(res_usage.ru_nivcsw);
        rv["rusage"]["minflt"] = Json::Int(res_usage.ru_minflt);
        rv["rusage"]["majflt"] = Json::Int(res_usage.ru_majflt);
    } else {
        rv["rusage"] = Json::Value(Json::nullValue);
    }
#endif

    return rv;
}

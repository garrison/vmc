#ifndef _VMC_RUN_INFORMATION_HPP
#define _VMC_RUN_INFORMATION_HPP

namespace RunInformation
{
    extern const char *compiler;

    namespace Precision
    {
        extern int digits;
        extern int min_exponent;
        extern int max_exponent;
    };

    extern const char *boost_version;
    extern const char *eigen_version;
};

#endif

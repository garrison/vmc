#ifndef _VMC_NOT_IMPLEMENTED_HPP
#define _VMC_NOT_IMPLEMENTED_HPP

#include <stdexcept>

/**
 * Exception for something that could be implemented but is currently not
 * supported by the code.
 *
 * This derives from std::invalid_argument so that it will result in a
 * ValueError in python code that is wrapped using Cython.
 */
class vmc_not_implemented : public std::invalid_argument
{
public:
    vmc_not_implemented (const char *why)
        : std::invalid_argument(why)
        {
        }
};

#endif

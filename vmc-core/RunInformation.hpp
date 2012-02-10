#ifndef _RUN_INFORMATION_HPP
#define _RUN_INFORMATION_HPP

#include <ctime>

#include <json/json.h>

class RunInformation
{
public:
    /**
     * This should be instantiated when the program begins.  Otherwise, the
     * "start_time" will be inaccurate.
     */
    RunInformation (void);

    /**
     * Return a json object with all information about this run.
     */
    Json::Value json_info (void) const;

private:
    std::time_t start_time;
};

#endif

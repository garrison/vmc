#ifndef _RUN_INFORMATION_HPP
#define _RUN_INFORMATION_HPP

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
};

#endif

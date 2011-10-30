#ifndef _SUBSYSTEM_HPP
#define _SUBSYSTEM_HPP

template<class T>
class Subsystem
// this is an abstract base class
{
public:
    virtual bool particle_is_within (const typename T::position_t &position) const = 0;

    virtual ~Subsystem (void)
	{
	}
};

#endif

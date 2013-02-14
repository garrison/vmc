from random import randint

cdef class RandomNumberGenerator(object):
    def __init__(self, name="boost::mt19937", seed=None):
        namebytes = name.encode("ascii")
        if not rng_name_is_valid(namebytes):
            raise ValueError("invalid RNG: {}".format(name))
        self.name = name
        if seed is None:
            self.seed = randint(0, <unsigned long>(-1))
        else:
            self.seed = seed
        self.autoptr = create_rng(namebytes, self.seed)

    def is_good(self):
        return self.autoptr.get() is not NULL

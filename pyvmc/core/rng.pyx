from random import randint

cdef class RandomNumberGenerator(object):
    def __init__(self):
        self.seed = randint(0, <unsigned long>(-1))
        self.autoptr = create_rng("boost::mt19937", self.seed)

    def is_good(self):
        return self.autoptr.get() is not NULL

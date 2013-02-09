from abc import ABCMeta
from collections import Hashable

# initial inspiration from http://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python
# I got a bit carried away, I think...

class ImmutableMetaclass(ABCMeta):
    def __init__(cls, name, bases, attrs):
        try:
            Immutable
        except NameError:
            # do nothing special while constructing the Immutable base class
            return
        if "__slots__" not in attrs:
            raise Exception("")
        if not isinstance(cls.__slots__, (list, tuple)):
            raise TypeError("")
        slots_set = frozenset(cls.__slots__)
        if len(cls.__slots__) != len(slots_set):
            raise Exception
        if "_immutable_slots" in attrs:
            if not isinstance(cls._immutable_slots, (list, tuple)):
                raise TypeError("")
            immutable_slots_set = frozenset(cls._immutable_slots)
            if len(cls._immutable_slots) != len(immutable_slots_set):
                raise Exception("")
            if not slots_set.issuperset(immutable_slots_set):
                raise Exception("")
        else:
            cls._immutable_slots = tuple(attrs["__slots__"])
            attrs["_immutable_slots"] = cls._immutable_slots
        super(ImmutableMetaclass, cls).__init__(name, bases, attrs)

class Immutable(Hashable):
    """Base class for various immutable objects.

    Implements various methods that an immutable object should have, including
    __eq__(), __ne__(), __hash__(), and __repr__().

    __slots__ must be defined, and must contain all attributes that the object
    might hold.

    _immutable_slots refers to the subset of __slots__ which are actually
    immutable.  If not given, it defaults to the value of __slots__.

    Attributes which are in __slots__ but not _immutable_slots should not
    affect the identity of the object; they are not considered in hashes or
    equality tests.  They exist so that if there is an expensive method of the
    object, its results can be cached in such an attribute.  You have to call
    object.__setattr__ to set any attribute in _immutable_slots, in order to
    get around the Immutable object's write protection.

    Never subclass the __init__() method.  Instead, init_validate() should be a
    function that has arguments that precisely match _immutable_slots.  Any
    assertions can be done here, and any default arguments can be given in the
    definition of init_validate().  init_validate() MUST return a tuple
    representing the intended (immutable) value of each attribute given in
    __slots__.  Each element of this tuple must be a hashable object.
    """

    __metaclass__ = ImmutableMetaclass

    __slots__ = ()

    def __init__(self, *args, **kwargs):
        args = self.init_validate(*args, **kwargs)
        assert isinstance(args, tuple)
        assert len(args) == len(self._immutable_slots)
        for slot, arg in zip(self._immutable_slots, args):
            assert isinstance(arg, Hashable)
            object.__setattr__(self, slot, arg)

    def __eq__(self, other):
        if self is other:
            return True
        return (self.__class__ == other.__class__ and
                all(getattr(self, attr) == getattr(other, attr)
                    for attr in self._immutable_slots))

    def __ne__(self, other):
        if self is other:
            return False
        return (self.__class__ != other.__class__ or
                any(getattr(self, attr) != getattr(other, attr)
                    for attr in self._immutable_slots))

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               ', '.join(repr(getattr(self, attr))
                                         for attr in self._immutable_slots))

    def __hash__(self):
        return hash(tuple(getattr(self, attr) for attr in self._immutable_slots))

    def __setattr__(self, name, value):
        raise TypeError("Immutable object")

    def __delattr__(self, name):
        raise TypeError("Immutable object")

    def init_validate(self, *args, **kwargs):
        while len(args) < len(self._immutable_slots):
            next_arg = self._immutable_slots[len(args)]
            try:
                args += (kwargs.pop(next_arg),)
            except KeyError:
                raise TypeError("argument not provided: {}".format(next_arg))
        if kwargs:
            raise TypeError("unexpected keyword argument: {}".format(kwargs.keys()[0]))
        if len(args) != len(self._immutable_slots):
            raise TypeError("Expected {0} arguments but received {1}".format(len(self._immutable_slots),
                                                                             len(args)))
        return args

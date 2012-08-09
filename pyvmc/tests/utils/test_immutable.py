from pyvmc.utils.immutable import Immutable

def test_immutable():
    class SomeImmutable(Immutable):
        __slots__ = ("one", "two")

    class SomeOtherImmutable(Immutable):
        __slots__ = ("one", "two", "_cached")
        _immutable_slots = ("one", "two")

    class DerivedImmutable(SomeImmutable):
        __slots__ = ("one", "two", "three")

    class TragicDerivedImmutable(SomeImmutable):
        __slots__ = ("one", "three")

    SomeImmutable(1, 2)
    SomeOtherImmutable(1, 2)
    DerivedImmutable(1, 2, 3)
    TragicDerivedImmutable(1, 3)

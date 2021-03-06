#include "Lattice.hpp"

bool LatticeSite::operator< (const LatticeSite &other) const
{
    if (this->n_dimensions() != other.n_dimensions())
        return bool(this->n_dimensions() < other.n_dimensions());
    if (this->basis_index != other.basis_index)
        return bool(this->basis_index < other.basis_index);

    unsigned int i = this->n_dimensions();
    for (;;) {
        if (i == 0)
            return false;
        --i;
        if (this->bs[i] != other.bs[i])
            return this->bs[i] < other.bs[i];
    }
}

Lattice::Lattice (const Lattice::DimensionVector &dimensions_, int basis_indices_)
    : dimensions(dimensions_),
      basis_indices(basis_indices_),
      m_total_sites(count_total_sites(dimensions_, basis_indices_)),
      stride(dimensions_.size())
{
    assert(dimensions.size() > 0);
    assert(basis_indices > 0);

    // set up strides
    unsigned int c = 1;
    for (unsigned int i = 0; i < dimensions.size(); ++i) {
        stride[i] = c;
        c *= dimensions[i];
    }
    basis_stride = c;

    // set up default move axes
    for (unsigned int i = 0; i < dimensions.size(); ++i) {
        MoveAxis m;
        m.bravais_site.resize(dimensions.size(), 0);
        m.bravais_site[i] = 1;
        m.basis_index = 0;
        move_axes.push_back(m);
    }
    if (basis_indices > 1) {
        MoveAxis m;
        m.bravais_site.resize(dimensions.size(), 0);
        m.basis_index = 1;
        move_axes.push_back(m);
    }
}

Lattice::Lattice (const std::initializer_list<int> &dimensions_, int basis_indices_)
    : Lattice(Lattice::DimensionVector(dimensions_), basis_indices_)
{
}

LatticeSite Lattice::operator[] (unsigned int n) const
{
    assert(n < total_sites());
    LatticeSite rv(this->n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        rv[i] = n % dimensions[i];
        n /= dimensions[i];
    }
    rv.basis_index = n;
#ifdef VMC_CAREFUL
    assert(site_is_valid(rv));
#endif
    return rv;
}

unsigned int Lattice::index (const LatticeSite &site) const
{
    assert(site_is_valid(site));

    unsigned int n = 0;
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        n += site[i] * stride[i];
    }
    n += site.basis_index * basis_stride;
#ifdef VMC_CAREFUL
    assert(site == this->index(n));
#endif
    return n;
}

bool Lattice::site_is_valid (const LatticeSite &site) const
{
    assert(site.n_dimensions() == n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        if (site[i] >= dimensions[i] || site[i] < 0)
            return false;
    }
    if (site.basis_index >= basis_indices || site.basis_index < 0)
        return false;
    return true;
}

void Lattice::asm_add_site_vector (LatticeSite &site, const BravaisSite &other) const
{
    assert(site.n_dimensions() == n_dimensions());
    assert(site.basis_index < basis_indices);
    assert(other.size() == n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i)
        site[i] += other[i];
}

void Lattice::asm_subtract_site_vector (LatticeSite &site, const BravaisSite &other) const
{
    assert(site.n_dimensions() == n_dimensions());
    assert(site.basis_index < basis_indices);
    assert(other.size() == n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i)
        site[i] -= other[i];
}

void Lattice::enforce_periodic_boundary (LatticeSite &site) const
{
    assert(site.n_dimensions() == n_dimensions());
    for (unsigned int dim = 0; dim < n_dimensions(); ++dim) {
        do_safe_modulus(site[dim], dimensions[dim]);
    }

    // this is often unnecessary ... should it be in a separate
    // function to be called before this one when needed?
    do_safe_modulus(site.basis_index, basis_indices);

    assert(site_is_valid(site));
}

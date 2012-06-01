#include "Lattice.hpp"

LatticeSite Lattice::site_from_index (unsigned int n) const
{
    BOOST_ASSERT(n < total_sites());
    LatticeSite rv(*this);
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        rv[i] = n % dimensions[i];
        n /= dimensions[i];
    }
    rv.basis_index = n;
#ifdef CAREFUL
    BOOST_ASSERT(site_is_valid(rv));
#endif
    return rv;
}

unsigned int Lattice::site_to_index (const LatticeSite &site) const
{
    BOOST_ASSERT(site_is_valid(site));

    unsigned int n = 0;
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        n += site[i] * offset[i];
    }
    n += site.basis_index * basis_offset;
#ifdef CAREFUL
    BOOST_ASSERT(site == site_from_index(n));
#endif
    return n;
}

bool Lattice::site_is_valid (const LatticeSite &site) const
{
    if (site.bravais_site().size() != n_dimensions())
        return false;
    for (unsigned int i = 0; i < n_dimensions(); ++i) {
        if (site[i] >= dimensions[i] || site[i] < 0)
            return false;
    }
    if (site.basis_index >= basis_indices || site.basis_index < 0)
        return false;
    return true;
}

phase_t Lattice::asm_add_site_vector (LatticeSite &site, const BravaisSite &other, const BoundaryConditions *bcs) const
{
    BOOST_ASSERT(site.bravais_site().size() == n_dimensions());
    BOOST_ASSERT(site.basis_index < basis_indices);
    BOOST_ASSERT(other.size() == n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i)
        site[i] += other[i];
    return enforce_boundary(site, bcs);
}

phase_t Lattice::asm_subtract_site_vector (LatticeSite &site, const BravaisSite &other, const BoundaryConditions *bcs) const
{
    BOOST_ASSERT(site.bravais_site().size() == n_dimensions());
    BOOST_ASSERT(site.basis_index < basis_indices);
    BOOST_ASSERT(other.size() == n_dimensions());
    for (unsigned int i = 0; i < n_dimensions(); ++i)
        site[i] -= other[i];
    return enforce_boundary(site, bcs);
}

phase_t Lattice::enforce_boundary (LatticeSite &site, const BoundaryConditions *bcs) const
{
    BOOST_ASSERT(site.bravais_site().size() == n_dimensions());
    BOOST_ASSERT(!bcs || bcs->size() == n_dimensions());
    phase_t phase_change = 1;
    for (unsigned int dim = 0; dim < n_dimensions(); ++dim) {
        while (site[dim] >= dimensions[dim]) {
            site[dim] -= dimensions[dim];
            if (bcs)
                phase_change *= (*bcs)[dim].phase();
        }
        while (site[dim] < 0) {
            site[dim] += dimensions[dim];
            if (bcs)
                phase_change /= (*bcs)[dim].phase();
        }
    }

    // this is often unnecessary ... should it be in a separate
    // function to be called before this one when needed?
    do_safe_modulus(site.basis_index, basis_indices);

    BOOST_ASSERT(site_is_valid(site));
    return phase_change;
}

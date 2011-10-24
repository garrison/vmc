#ifndef _SWAPPED_SYSTEM_HPP
#define _SWAPPED_SYSTEM_HPP

// fixme: this is polluting the namespace.  isn't there something like this in <algorithm> already, anyway?
template <typename T>
int vector_find (const std::vector<T> &vec, const T &val)
{
    for (unsigned int i = 0; i < vec.size(); ++i) {
	if (vec[i] == val)
	    return i;
    }
    return -1;
}

template <class T>
class SwappedSystem
{
private:
    CeperlyMatrix<typename T::amplitude_t> phibeta1, phibeta2;
    const Subsystem<T> *subsystem;
    std::vector<unsigned int> copy1_subsystem_indices, copy2_subsystem_indices;
    bool phibeta1_is_dirty, phibeta2_is_dirty;
    bool initialized;
public:
    SwappedSystem (const Subsystem<T> &subsystem_)
	: subsystem(&subsystem_),
	  phibeta1_is_dirty(false),
	  phibeta2_is_dirty(false),
	  initialized(false)
	{
	}

    void initialize (const typename T::Arguments &r1, const typename T::Arguments &r2,
		     const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2)
	{
	    BOOST_ASSERT(!initialized);

	    unsigned int N = r1.size();
	    for (unsigned int i = 0; i < N; ++i) {
		if (subsystem->particle_is_within(r1[i]))
		    copy1_subsystem_indices.push_back(i);
		if (subsystem->particle_is_within(r2[i]))
		    copy2_subsystem_indices.push_back(i);
	    }

	    Eigen::Matrix<typename T::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> phibeta1_tmp(phialpha1.get_matrix()), phibeta2_tmp(phialpha2.get_matrix());
	    unsigned int particles_to_swap = std::min(copy1_subsystem_indices.size(), copy2_subsystem_indices.size());
	    std::cerr << "swapping " << particles_to_swap << " particles" << std::endl;
	    for (unsigned int i = 0; i < particles_to_swap; ++i) {
		phibeta1_tmp.row(copy1_subsystem_indices[i]).swap(phibeta2_tmp.row(copy2_subsystem_indices[i]));
	    }
	    phibeta1 = phibeta1_tmp;
	    phibeta2 = phibeta2_tmp;

	    initialized = true;

#if 1
	    verify_phibetas(r1, r2, phialpha1, phialpha2);
#endif
	}

    void verify_phibetas (const typename T::Arguments &r1, const typename T::Arguments &r2,
			  const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2)
	{
	    BOOST_ASSERT(initialized);

	    unsigned int N = r1.size();

	    // verify that the system index arrays have everything they need (and no duplicates!)
	    unsigned int c1 = 0, c2 = 0;
	    for (unsigned int i = 0; i < N; ++i) {
		bool b1 = vector_find(copy1_subsystem_indices, i) != -1;
		bool b2 = vector_find(copy2_subsystem_indices, i) != -1;
		if (b1) ++c1;
		if (b2) ++c2;
		BOOST_ASSERT(b1 == subsystem->particle_is_within(r1[i]));
		BOOST_ASSERT(b2 == subsystem->particle_is_within(r2[i]));
	    }
	    BOOST_ASSERT(c1 == copy1_subsystem_indices.size());
	    BOOST_ASSERT(c2 == copy2_subsystem_indices.size());

	    // fixme: these next few lines are copy-pasted from above ...
	    Eigen::Matrix<typename T::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> phibeta1_tmp(phialpha1.get_matrix()), phibeta2_tmp(phialpha2.get_matrix());
	    unsigned int particles_to_swap = std::min(copy1_subsystem_indices.size(), copy2_subsystem_indices.size());
	    for (unsigned int i = 0; i < particles_to_swap; ++i) {
		phibeta1_tmp.row(copy1_subsystem_indices[i]).swap(phibeta2_tmp.row(copy2_subsystem_indices[i]));
	    }

#if 1
	    for (unsigned int i = 0; i < N; ++i) {
		    if (phibeta1_tmp(i, 0) != phibeta1.get_matrix()(i, 0))
			std::cerr << "Alert1! " << i << std::endl;
		    if (phibeta2_tmp(i, 0) != phibeta2.get_matrix()(i, 0))
			std::cerr << "Alert2! " << i << std::endl;
	    }
#endif
	    for (unsigned int i = 0; i < N; ++i) {
		for (unsigned int j = 0; j < N; ++j) {
		    BOOST_ASSERT(phibeta1_tmp(i, j) == phibeta1.get_matrix()(i, j));
		    BOOST_ASSERT(phibeta2_tmp(i, j) == phibeta2.get_matrix()(i, j));
		}
	    }
	}

    void update (int copy, unsigned int index,
		 const typename T::Arguments &r1, const typename T::Arguments &r2,
		 const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2,
		 typename T::amplitude_t &phibeta1_ratio_multiplier, typename T::amplitude_t &phibeta2_ratio_multiplier)
	{
	    BOOST_ASSERT(initialized);
	    BOOST_ASSERT(copy == 1 || copy == 2);

	    const typename T::Arguments &r_this = (copy == 1) ? r1 : r2;
	    const CeperlyMatrix<typename T::amplitude_t> &phialpha_this = (copy == 1) ? phialpha1 : phialpha2;
	    const CeperlyMatrix<typename T::amplitude_t> &phialpha_that = (copy == 1) ? phialpha2 : phialpha1;
	    CeperlyMatrix<typename T::amplitude_t> &phibeta_this = (copy == 1) ? phibeta1 : phibeta2;
	    CeperlyMatrix<typename T::amplitude_t> &phibeta_that = (copy == 1) ? phibeta2 : phibeta1;
	    std::vector<unsigned int> &subsystem_indices_this = (copy == 1) ? copy1_subsystem_indices : copy2_subsystem_indices;
	    std::vector<unsigned int> &subsystem_indices_that = (copy == 1) ? copy2_subsystem_indices : copy1_subsystem_indices;
	    typename T::amplitude_t &phibeta_this_ratio_multiplier = (copy == 1) ? phibeta1_ratio_multiplier : phibeta2_ratio_multiplier;
	    typename T::amplitude_t &phibeta_that_ratio_multiplier = (copy == 1) ? phibeta2_ratio_multiplier : phibeta1_ratio_multiplier;
	    bool &phibeta_this_is_dirty = (copy == 1) ? phibeta1_is_dirty : phibeta2_is_dirty;
	    bool &phibeta_that_is_dirty = (copy == 1) ? phibeta2_is_dirty : phibeta1_is_dirty;

	    int si = vector_find(subsystem_indices_this, index);
	    bool is_in_subsystem_now = subsystem->particle_is_within(r_this[index]);
	    if (si == -1) {
		// particle was not in subsystem before
		if (is_in_subsystem_now) {
		    // particle entered subsystem
		    subsystem_indices_this.push_back(index);
		    if (subsystem_indices_this.size() <= subsystem_indices_that.size()) {
			std::cerr << "A";
			// we can (and should) actually perform a swap
			si = subsystem_indices_this.size() - 1;
			phibeta_this.update_row(subsystem_indices_this[si], phialpha_that.get_matrix().row(subsystem_indices_that[si]));
			phibeta_that.update_row(subsystem_indices_that[si], phialpha_this.get_matrix().row(subsystem_indices_this[si]));
			phibeta_that_ratio_multiplier *= phibeta_that.calculate_determinant_ratio();
			phibeta_that_is_dirty = true;
		    } else {
			std::cerr << "B";
			// we cannot perform a swap due to subsystem count mismatch
			phibeta_this.update_row(index, phialpha_this.get_matrix().row(index));
		    }
		} else {
		    std::cerr << "-";
		    // still not in subsystem.  trivial update.
		    phibeta_this.update_row(index, phialpha_this.get_matrix().row(index));
		}
		phibeta_this_ratio_multiplier *= phibeta_this.calculate_determinant_ratio();
		phibeta_this_is_dirty = true;
	    } else {
		// particle was in subsystem before
		if (!is_in_subsystem_now) {
		    // particle left subsystem
		    if ((unsigned int) si < subsystem_indices_that.size()) {
			if (subsystem_indices_this.size() <= subsystem_indices_that.size()) {
			    // we will have to leave this particle unpaired
			    std::cerr << "D";
			    // we need to unswap
			    phibeta_this.update_row(index, phialpha_this.get_matrix().row(index));
			    phibeta_that.update_row(subsystem_indices_that[si], phialpha_that.get_matrix().row(subsystem_indices_that[si]));
			    phibeta_that_ratio_multiplier *= phibeta_that.calculate_determinant_ratio();
			    phibeta_that_is_dirty = true;
			    // and now we rearrange the subsystem indices so that they are paired correctly
			    size_t last_common_index = std::min(subsystem_indices_this.size(), subsystem_indices_that.size()) - 1;
//			    std::cerr << " ... " << si << " " << last_common_index << "   " << subsystem_indices_this[si] << " (" << index << ") " << subsystem_indices_this[last_common_index] << "   " << subsystem_indices_that[si] << " " << subsystem_indices_that[last_common_index] << "     " << subsystem_indices_this.size() - 1 << " " << subsystem_indices_this[subsystem_indices_this.size() - 1] << "        " << subsystem_indices_this.size() << " " << subsystem_indices_that.size() << std::endl;
			    std::swap(subsystem_indices_that[si], subsystem_indices_that[last_common_index]);
			    subsystem_indices_this[si] = subsystem_indices_this[last_common_index];
			    subsystem_indices_this[last_common_index] = subsystem_indices_this[subsystem_indices_this.size() - 1];
			    subsystem_indices_this.pop_back();
			} else {
			    // there's no reason to leave this particle unpaired.  pair it!
			    std::cerr << "G";
			    subsystem_indices_this[si] = subsystem_indices_this[subsystem_indices_this.size() - 1];
			    subsystem_indices_this.pop_back();
			    phibeta_this.update_row(index, phialpha_this.get_matrix().row(index)); // this duplicates above
			    phibeta_this_ratio_multiplier *= phibeta_this.calculate_determinant_ratio();
			    phibeta_this.finish_row_update();
			    phibeta_this.update_row(subsystem_indices_this[si], phialpha_that.get_matrix().row(subsystem_indices_that[si]));
			    phibeta_that.update_row(subsystem_indices_that[si], phialpha_this.get_matrix().row(subsystem_indices_this[si]));
			    phibeta_that_ratio_multiplier *= phibeta_that.calculate_determinant_ratio();
			    phibeta_that_is_dirty = true;
			}
		    } else {
			std::cerr << "E";
			// this particle had never been swapped to begin with, due to count mismatch
			phibeta_this.update_row(index, phialpha_this.get_matrix().row(index));
			// rearrange it from subsystem indices
			subsystem_indices_this[si] = subsystem_indices_this[subsystem_indices_this.size() - 1];
			subsystem_indices_this.pop_back();
		    }
		    phibeta_this_ratio_multiplier *= phibeta_this.calculate_determinant_ratio();
		    phibeta_this_is_dirty = true;
		} else {
		    // particle stayed in subsystem.  mostly trivial update.
		    if ((unsigned int) si < subsystem_indices_that.size()) {
			std::cerr << "_";
			phibeta_that.update_row(subsystem_indices_that[si], phialpha_this.get_matrix().row(index));
			phibeta_that_ratio_multiplier *= phibeta_that.calculate_determinant_ratio();
			phibeta_that_is_dirty = true;
		    } else {
			std::cerr << "^";
			phibeta_this.update_row(index, phialpha_this.get_matrix().row(index));
			phibeta_this_ratio_multiplier *= phibeta_this.calculate_determinant_ratio();
			phibeta_this_is_dirty = true;
		    }
		}
	    }
#if 1
	    verify_phibetas(r1, r2, phialpha1, phialpha2); // fixme: enable this only for debugging!!
#endif
	}

    void finish_update (void)
	{
	    BOOST_ASSERT(initialized);
#ifdef DEBUG
	    std::cerr << "finished swap-system update" << std::endl;
#endif

	    if (phibeta1_is_dirty)
		phibeta1.finish_row_update();
	    if (phibeta2_is_dirty)
		phibeta2.finish_row_update();
	    phibeta1_is_dirty = false;
	    phibeta2_is_dirty = false;
	}

    int calculate_subsystem_particle_change (int copy, unsigned int index, const typename T::Arguments &r)
	{
	    BOOST_ASSERT(initialized);
	    BOOST_ASSERT(copy == 1 || copy == 2);

	    const std::vector<unsigned int> &subsystem_indices = (copy == 1) ? copy1_subsystem_indices : copy2_subsystem_indices;

	    int rv = 0;
	    if (subsystem->particle_is_within(r[index]))
		rv += 1;
	    if (vector_find(subsystem_indices, index) != -1)
		rv -= 1;
	    return rv;
	}

    const CeperlyMatrix<typename T::amplitude_t> & get_phibeta1 (void) const
	{
	    BOOST_ASSERT(initialized);
	    return phibeta1;
	}

    const CeperlyMatrix<typename T::amplitude_t> & get_phibeta2 (void) const
	{
	    BOOST_ASSERT(initialized);
	    return phibeta2;
	}

    unsigned int get_N_subsystem1 (void) const
	{
	    return copy1_subsystem_indices.size();
	}

    unsigned int get_N_subsystem2 (void) const
	{
	    return copy2_subsystem_indices.size();
	}

private:
    // disable default constructor
    SwappedSystem (void);
};

#endif

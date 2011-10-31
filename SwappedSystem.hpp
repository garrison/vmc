#ifndef _SWAPPED_SYSTEM_HPP
#define _SWAPPED_SYSTEM_HPP

#include <vector>
#include <list>

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

// fixme: also polluting namespace
template <typename T>
void no_duplicate_append (std::list<T> &p, const T &val)
{
    for (typename std::list<T>::const_iterator i = p.begin(); i != p.end(); ++i) {
	if (*i == val)
	    return;
    }
    p.push_back(val);
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
    SwappedSystem (const Subsystem<T> *subsystem_)
	: subsystem(subsystem_),
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

	    if (copy1_subsystem_indices.size() == copy2_subsystem_indices.size())
		reinitialize_phibetas(r1, r2, phialpha1, phialpha2);

	    initialized = true;
	}

private:
    void reinitialize_phibetas (const typename T::Arguments &r1, const typename T::Arguments &r2,
				const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2)
	{
	    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());

	    Eigen::Matrix<typename T::amplitude_t, Eigen::Dynamic, Eigen::Dynamic> phibeta1_tmp(phialpha1.get_matrix()), phibeta2_tmp(phialpha2.get_matrix());
	    const unsigned int particles_to_swap = copy1_subsystem_indices.size();
#ifdef DEBUG
	    std::cerr << "swapping " << particles_to_swap << " particles" << std::endl;
#endif
	    for (unsigned int i = 0; i < particles_to_swap; ++i) {
		phibeta1_tmp.row(copy1_subsystem_indices[i]).swap(phibeta2_tmp.row(copy2_subsystem_indices[i]));
	    }
	    phibeta1 = phibeta1_tmp;
	    phibeta2 = phibeta2_tmp;

#if 1
	    verify_phibetas(r1, r2, phialpha1, phialpha2);
#endif
	    // variables needed only for verification (unused otherwise)
	    (void) r1;
	    (void) r2;
	}

public:
    void verify_phibetas (const typename T::Arguments &r1, const typename T::Arguments &r2,
			  const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2) const
	{
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

#ifdef CAREFUL
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

    void update (int index1, int index2,
		 const typename T::Arguments &r1, const typename T::Arguments &r2,
		 const CeperlyMatrix<typename T::amplitude_t> &phialpha1, const CeperlyMatrix<typename T::amplitude_t> &phialpha2)
	{
	    BOOST_ASSERT(initialized);
	    BOOST_ASSERT(index1 >= -1 && index1 < (int) r1.size());
	    BOOST_ASSERT(index2 >= -1 && index2 < (int) r1.size());

	    const bool phibetas_were_up_to_date = (get_N_subsystem1() == get_N_subsystem2());

	    std::list<unsigned int> touched_1, touched_2;

	    if (index1 != -1)
		update_subsystem_indices(subsystem, index1, r1, touched_1, touched_2, copy1_subsystem_indices, copy2_subsystem_indices);
	    if (index2 != -1)
		update_subsystem_indices(subsystem, index2, r2, touched_2, touched_1, copy2_subsystem_indices, copy1_subsystem_indices);

	    if (get_N_subsystem1() != get_N_subsystem2()) {
		// we can't perform a swap, so we don't even attempt to update
		// the phibeta's
#ifdef DEBUG_SWAPPEDSYSTEM
		std::cerr << "-";
#endif
		return;
	    }

	    if (phibetas_were_up_to_date && touched_1.size() < 2) { // FIXME!!
		// use CeperlyMatrix update steps
#ifdef DEBUG_SWAPPEDSYSTEM
		std::cerr << ' ' << touched_1.size() << touched_2.size();
#endif

		update_touched_rows(subsystem, r1, touched_1, copy1_subsystem_indices, copy2_subsystem_indices, phialpha1, phibeta1, phibeta2);
		update_touched_rows(subsystem, r2, touched_2, copy2_subsystem_indices, copy1_subsystem_indices, phialpha2, phibeta2, phibeta1);
	    } else {
		// reinitialize CeperlyMatrices
#ifdef DEBUG_SWAPPEDSYSTEM
		std::cerr << "*";
#endif
		reinitialize_phibetas(r1, r2, phialpha1, phialpha2);
	    }

#if 1
	    verify_phibetas(r1, r2, phialpha1, phialpha2); // fixme: enable this only for debugging!!
#endif
	}

    // fixme: finish_update() function ??  will save some time on swapa,sign calculation

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
	    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());
	    return phibeta1;
	}

    const CeperlyMatrix<typename T::amplitude_t> & get_phibeta2 (void) const
	{
	    BOOST_ASSERT(initialized);
	    BOOST_ASSERT(copy1_subsystem_indices.size() == copy2_subsystem_indices.size());
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

    static void update_subsystem_indices (const Subsystem<T> *subsystem, unsigned int index, const typename T::Arguments &r,
					  std::list<unsigned int> &touched_this, std::list<unsigned int> &touched_that,
					  std::vector<unsigned int> &subsystem_indices_this, std::vector<unsigned int> &subsystem_indices_that)
	{
	    int si = vector_find(subsystem_indices_this, index);
	    no_duplicate_append(touched_this, index);
	    if (subsystem->particle_is_within(r[index])) {
		if (si == -1) {
		    // particle entered subsystem
#if 0
		    si = (int) subsystem_indices_this.size();
#endif
		    subsystem_indices_this.push_back(index);
		}
	    } else {
		if (si != -1) {
		    // particle left subsystem
		    subsystem_indices_this[si] = subsystem_indices_this[subsystem_indices_this.size() - 1];
		    no_duplicate_append(touched_this, subsystem_indices_this[si]);
		    if (si < (int) subsystem_indices_that.size())
			no_duplicate_append(touched_that, subsystem_indices_that[si]);
		    subsystem_indices_this.pop_back();
#if 0
		    si = -1;
#endif
		}
	    }
	}

    static void update_touched_rows (const Subsystem<T> *subsystem, const typename T::Arguments &r,
				     const std::list<unsigned int> &touched_this,
				     const std::vector<unsigned int> &subsystem_indices_this,
				     const std::vector<unsigned int> &subsystem_indices_that,
				     const CeperlyMatrix<typename T::amplitude_t> &phialpha_this,
				     CeperlyMatrix<typename T::amplitude_t> &phibeta_this,
				     CeperlyMatrix<typename T::amplitude_t> &phibeta_that)
	{
	    for (std::list<unsigned int>::const_iterator i = touched_this.begin(); i != touched_this.end(); ++i) {
		if (subsystem->particle_is_within(r[*i])) {
		    // it goes in the other phibeta
		    const int pair_index = vector_find(subsystem_indices_this, *i);
		    BOOST_ASSERT(pair_index != -1);
		    phibeta_that.update_row(subsystem_indices_that[pair_index], phialpha_this.get_matrix().row(*i));
		    phibeta_that.finish_row_update();
		} else {
		    // it goes in this phibeta
		    phibeta_this.update_row(*i, phialpha_this.get_matrix().row(*i));
		    phibeta_this.finish_row_update();
		}
	    }
	}

};

#endif

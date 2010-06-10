#ifndef OB_PERMUTATION_H
#define OB_PERMUTATION_H

#include <openbabel/babelconfig.h>

#include <map>
#include <vector>
#include <algorithm>

#include <Eigen/Core>

namespace OpenBabel {

  class OBMol;

  struct OBAPI Permutation
  {
    /**
     * Default constructor.
     */
    Permutation()
    {}

    /**
     * Constructor taking a @p map as argument.
     */
    Permutation(const std::vector<unsigned int> &_map) : map(_map)
    {}

    /**
     * Constructor taking a @p matrix as argument.
     */
    Permutation(const Eigen::MatrixXi &matrix)
    {
      setMatrix(matrix);
    }

    /**
     * Copy constructor.
     */
    Permutation(const Permutation &other)
    {
      this->map = other.map;
    }

    /**
     * Print the permutation to std::cout in shortened notation.
     */
    void print() const
    {
      std::vector<unsigned int>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i)
        std::cout << *i << " ";
      std::cout << std::endl;    
    }

    Permutation apply(const Permutation &input) const
    {
      Permutation p;
      if (input.map.size() != map.size())
        return p;
      p.map.resize(map.size());

      unsigned int index = 0;
      std::vector<unsigned int>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i, ++index) {
        p.map[index] = input.map.at(*i-1);
      }

      return p;   
    }

    /**
     * Compute the permutation matrix for this Permutation.
     */
    Eigen::MatrixXi matrix() const
    {
      int n = map.size();
      Eigen::MatrixXi P = Eigen::MatrixXi::Zero(n, n);

      // make sure we can create matrices for contracted permutations 
      // (i.e. permutations where the values don't go from 1 to N, there
      // may be gaps)
      std::map<unsigned int, unsigned int> renum;
      unsigned int j = 0;
      unsigned int k = 1;
      while (renum.size() < n) {
        if (std::find(map.begin(), map.end(), k) != map.end()) {
          renum[k] = j;
          j++;
        }
      
        k++;
      }

      for (int i = 0; i < n; ++i) {
        //P(i, renum[map.at(i)]) = 1;
        P(renum[map.at(i)], i) = 1;
      }
      /*
      for (int i = 0; i < n; ++i) {
        //P(map.at(i)-1, i) = 1;
        P(i, map.at(i)-1) = 1;
      }
       */


      return P;
    }

    /**
     * Set the permutation matrix for this permutation. The matrix is
     * converted to a shortened notation permutation.
     */
    void setMatrix(const Eigen::MatrixXi &m)
    {
      if (m.rows() != m.cols())
        return;

      int size = m.rows();
      map.resize(size);

      for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
          if (m(i,j))
            map[i] = j + 1;
        }
    }

    unsigned long NumInversions() const
    {
      std::vector<unsigned int> invVec; // the inversion vector
      std::vector<unsigned int>::const_iterator i, j;
      for (i = map.begin(); i != map.end(); ++i) {
        int e = 0; // ith element
        // loop over elements to the right
        for (j = i; j != map.end(); ++j)
          // increment e if element to the right is lower
          if (*j < *i)
            e++;

        invVec.push_back(e);
      }

      unsigned long sum = 0;
      for (std::vector<unsigned int>::iterator k = invVec.begin(); k != invVec.end(); ++k)
        sum += *k;

      return sum;
    }


    /**
     * Multiply two permutation matrices and return the resulting permutation.
     * This is the same as applying the two permutations consecutively.
     */
    Permutation operator*(const Permutation &rhs) const
    {
      return apply(rhs);
    }

    /**
     * Equality operator. Two permutations are equal if their shortened notations 
     * are the same (compared element-by-element).
     */
    bool operator==(const Permutation &rhs) const
    {
      if (map.size() != rhs.map.size())
        return false;

      std::vector<unsigned int>::const_iterator i1 = map.begin();
      std::vector<unsigned int>::const_iterator i2 = rhs.map.begin();
      for (; i1 != map.end(); ++i1, ++i2)
        if (*i1 != *i2)
          return false;

      return true;  
    }

    std::vector<unsigned int> map; //!< The actual mapping in shortened notation
  };

  struct OBAPI PermutationGroup
  {
    /**
     * Default constructor.
     */
    PermutationGroup()
    {}

    /**
     * Constructor taking vector of permutations as argument.
     */
    PermutationGroup(const std::vector<Permutation> &_permutations) : permutations(_permutations)
    {}

    /**
     * @return The number of permutations in this group.
     */
    unsigned int size() const
    {
      return permutations.size();
    }

    /**
     * Add permutation @p p to this permutation group.
     */
    void add(const Permutation &p)
    {
      permutations.push_back(p);
      //invIndexes.push_back(p.NumInversions());
    }

    /**
     * @return A constant reference to the @p index-th permutation.
     */
    const Permutation& at(unsigned int index) const
    {
      return permutations.at(index);
    }

    /**
     * @return True if this permutation group contains permutation @p p.
     */
    bool contains(const Permutation &p) const
    {
      //if (std::find(invIndexes.begin(), invIndexes.end(), p.NumInversions()) != invIndexes.end())
      //  return true;
      std::vector<Permutation>::const_iterator i;
      for (i = permutations.begin(); i != permutations.end(); ++i)
        if (*i == p)
          return true;
      return false; 
    }

    std::vector<Permutation> permutations; //!< The actual permutations in the group
    //std::vector<unsigned long> invIndexes; //!< The inversion indexes
  };
  
  OBAPI PermutationGroup findAutomorphisms(OpenBabel::OBMol *obmol, const std::vector<unsigned int> &symClasses);

} // namespace

#endif

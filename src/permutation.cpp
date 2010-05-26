#include <openbabel/permutation.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>

#include <openbabel/bliss/graph.hh>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

namespace OpenBabel {

  void callback(void *param, unsigned int n, const unsigned int *aut)
  {
    // construct vector with elements
    std::vector<unsigned int> elements;
    //cout << "Generator: ";
    for (unsigned int i = 0; i < n; ++i) {
      //cout << aut[i] + 1 << " ";
      elements.push_back(aut[i]+1);
    }
    //cout << endl;

    // add the generator
    PermutationGroup *generators = static_cast<PermutationGroup*>(param);
    generators->add(Permutation(elements));
  }

  /**
   * Add the permutations' inverses to group @p G. These inverses automatically
   * belong to the automorphism group without the need of further testing.
   */
  void addInverses(PermutationGroup &G)
  {
    //cout << "addInverses..." << endl;
    unsigned int size = G.permutations.size();
    for (unsigned int i = 0; i < size; ++i) {
      int n = G.at(i).map.size();
      Eigen::MatrixXi P = Eigen::MatrixXi::Zero(n, n);
      for (int j = 0; j < n; ++j) {
        P(j, G.at(i).map.at(j)-1) = 1;
      }
 
      Permutation inv_p(P.transpose());
      //cout << "matrix:" << endl; cout << G.permutations.at(i).matrix() << endl;
      //inv_p.print();
      if (!G.contains(inv_p)) {
        //cout << "--------------> found inverse" << endl;
        G.add(inv_p); 
      }
    }
  }

  /**
   * Add the all permutation products to group @p G. These inverses automatically
   * belong to the automorphism group without the need of further testing.
   */
  void addProducts(PermutationGroup &G)
  {
    //cout << "addProducts..." << endl;
    for (unsigned int i = 0; i < G.permutations.size(); ++i) {
      for (unsigned int j = 0; j < G.permutations.size(); ++j) {
        if (i >= j)
          continue;

        Permutation p = G.permutations.at(i) * G.permutations.at(j);
        //p.print();
        if (!G.contains(p)) {
          //cout << "-------------> found product" << endl;
          G.add(p);
        }
      }
    }
  }

  PermutationGroup findAutomorphisms(OpenBabel::OBMol *obmol, const std::vector<unsigned int> &symClasses)
  {
    // construct the bliss graph
    bliss::Graph g;
    std::map<OpenBabel::OBAtom*, unsigned int> atom2index;
    FOR_ATOMS_OF_MOL (atom, obmol) {
      atom2index[&*atom] = g.add_vertex(symClasses.at(atom->GetIndex()));
    }
    FOR_BONDS_OF_MOL (bond, obmol) {
      g.add_edge(atom2index[bond->GetBeginAtom()], atom2index[bond->GetEndAtom()]);
    }

    // use bliss to get the automorphism group generators
    PermutationGroup generators;
    bliss::Stats stats;
    g.find_automorphisms(stats, &callback, &generators);

    //cout << "# Automorphisms:" << stats.group_size_approx << endl;
    unsigned long nAut = stats.group_size_approx;

    // construct the automorphism group
    PermutationGroup G;

    // add identity permutation
    std::vector<unsigned int> eElements;
    for (unsigned int i = 0; i < obmol->NumAtoms(); ++i)
      eElements.push_back(i+1);
    Permutation e(eElements);
    if (!G.contains(e))
      G.add(e);

    // add the generators
    for (unsigned int i = 0; i < generators.size(); ++i) {
      if (!G.contains(generators.at(i)))
        G.add(generators.at(i));
    }

    // loop and add inverses and products until no more automorphisms are found
    unsigned int counter = 0;
    unsigned int lastSize = 0;
    while (G.permutations.size() != lastSize) {
      lastSize = G.permutations.size();

      addInverses(G);
      addProducts(G);
    
      counter++;
      if (counter > 1000)
        break;
    }

    /*
    if (nAut != G.size())
      cout << "ERROR: Not all " << nAut << " automorphisms are found!" << endl;
    else 
      cout << "SUCCESS: all " << nAut << " automorphisms are found!" << endl;
      */
    
    return G;
  }

} // namespace

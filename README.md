# Computational Simulation of Multi-Body Interactions with O(N) Scaling

N-body employs a sequential multipole algorithm designed to accurately calculate the interparticle gravitational potential (1/r) for a system of N particles, achieving precision within round-off error limits. Notably, its computational efficiency scales linearly with the particle count (N), denoted as O(N). This scalability facilitates simulations involving a substantial number of atoms, surpassing 500,000 or even more. Realistic dynamic simulations often demand thousands to millions of time steps. The algorithm's linear computational cost for potential calculations enables a significantly increased number of time steps within a shorter duration compared to previous methods.

The multipole-expansion algorithm strategically calculates the potential by dividing it into two components: the contributions from nearby particles, referred to as 'near-field potentials,' and the potential arising from particles at a distance, known as 'far-field potentials.' Direct evaluation is employed to compute the near-field potential, while the far-field potential is expressed as an absolutely convergent Taylor series. To achieve a specified precision (e), the number of terms (p) in the Taylor series can be predetermined, ensuring a particular level of accuracy in the computation.

The objective of the multipole computation is to efficiently compute the far-field with a computational effort of O(N). This goal is realized through a two-phase process, involving the propagation of values both upwards and downwards within the computational tree.

During the first (upward) phase, the far-field interactions are calculated by shifting the multipole expansions of child nodes to the center of the parent node and summing the results. This computation spans all cubes at all levels, starting from the finest level of refinement (the leaf nodes). The multipole expansion of the field attributed to each particle is determined concerning the center of each leaf node and aggregated. Progressing upwards, the expansions of the eight child nodes are shifted relative to the parent's center and summed. This computation is executed at each level, propagating values upward through the tree.

By the end of this phase, each node at every level possesses a multipole-expansion representation of the field arising from all particles within its boundaries. Importantly, these expansions remain valid only outside the region defined by the near field of each cube.

In the second (downward) phase, we generate local expansions for all cubes at every level, commencing from the coarsest level (root node). The far-field multipole expansions of all nodes within the interactive field of the node under consideration are transformed into local expansions relative to the node's center and aggregated. Subsequently, the expansion of the parent cube is shifted to the center of the current cube under consideration and combined with the local expansions derived from the interactive field.

It's important to note that the root node lacks a parent or interactive field, making its local expansion defined as zero. The local expansions are then systematically propagated down the tree. These expansions encapsulate the field contributions from all particles in the system that are not contained within the current node or its closest neighbors.

Ultimately, for each cube i at the finest level n, we assess the local expansion for every particle within cube i. Each local expansion is characterized by the coefficients of a p-term Taylor series. The potential is obtained by directly evaluating this series at each particle. By summing the local expansion at each leaf node, we derive the potential or force arising from all particles beyond the nearest neighbor cubes.

The interactions of each particle in cube i are directly computed and aggregated to produce the potential attributed to the nearest neighbors. Finally, combining the far-field and near-field potentials through summation yields the desired potential.

A detailed discussion and a formal description of the multipole-expansion algorithm can be found in the following references:

  "Computational Simulation of Multi-Body Interactions with O(N) Scaling",
  Cristy and David A. Pensak, Central Research and Development,
  E.I. Du Pont de Nemours, Technical Report, June 1991.

  "An O(N) Algorithm for Three-dimensional N-body Simulations", Feng Zhao,
  MIT Artificial Intelligence Lab, Technical Report 995, 1987.

The formula for the lexicographical ordering of a tetrahedral array was
derived by Richard Merrifield.

## Build and run the multipole program

To build and run the multipole program, type

```
  ./configure
  make
  multipole
```

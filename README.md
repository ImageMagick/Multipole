# Computational Simulation of Multi-Body Interactions with O(N) Scaling

N-body implements a sequential multipole algorithm that solves the
interparticle gravitational potential 1/r of N particles to within
round-off error and whose computational effort grows only linearly with
increasing N, i.e. O(N).  This makes possible simulations of a large
number of atoms on the order of 500,000 or more.  Typically, thousands
to millions of time steps are required for a realistic dynamic
simulation.  The linear nature of the cost of computing the potential
with the sequential multipole algorithm allows many more time-steps, in
less time, than previously possible.

The multipole-expansion algorithm computes the potential by
partitioning them into two parts: close-by particles (`near-field
potentials') and the potential of far-away particles (`far-field
potentials').  The near-field potential is computed using direct
evaluation.  The far-field potential is represented as an absolutely
convergent Taylor series.  Given a required precision e, the number of
terms p in the Taylor series can be determined in advance to obtain a
particular accuracy.

The goal of the multipole computation is to compute the far-field with
O(N) work.  This is achieved by propagating values, first upwards
through the computational tree, then downward in two distinct phases.
In the first (upward) phase the far-field interactions are computed by
shifting the multipole expansions of the child nodes to the center of
the parent node and adding the results.  The far-field interaction is
calculated for all cubes at all levels, beginning at the finest level
of refinement (the leaf nodes).  The multipole-expansion of the field
due to each particle is calculated relative to the center of each leaf
node and summed.  Proceeding upward, the expansions of the eight child
nodes are shifted relative to the parent's center and summed.  This
computation is performed at each level propagating values upward
through the tree.  At the end of this phase of the algorithm, each node
at every level has a multipole-expansion representation of the field
due to all particles within its boundaries.  These expansions are only
valid outside the region defined by the near field of each cube.

In the second (downward) phase, we form the local expansions for all
cubes at all levels beginning at the coarsest level (root node).  The
far-field multipole-expansions of all nodes in the interactive field of
the node under consideration are converted into local expansions
relative to the node's center and are summed.  Next, shift the
expansion of the parent cube to the center of the current cube under
consideration and sum with the local expansions computed from the
interactive field.  Note, the root node does not have a parent or
interactive field so it's local expansion is defined to be zero.  The
local expansions are propagated down the tree.  They describe the field
due to all particles in the system that are not contained in the
current node or its nearest neighbors.

Finally, for each cube i at the finest level n, we evaluate the local
expansion for each particle in cube i.  Each local expansion is
described by the coefficients of a p-term Taylor series.  Direct
evaluation of this series at a particle yields the potential.
Summating the local expansion at each leaf node generates the potential
or force due to all particles outside the nearest neighbor cubes. The
interactions of each particle in cube i are computed directly and
summed to generate the potential due to the nearest neightbors.
Finally summating far-field and near-field potential gives us the
desired potential.

A discussion and formal description of the multipole-expansion algorithm
can be found in these references:

  "Computational Simulation of Multi-Body Interactions with O(N) Scaling",
  Cristy and David A. Pensak, Central Research and Development,
  E.I. Du Pont de Nemours, Technical Report, June 1991.

  "An O(N) Algorithm for Three-dimensional N-body Simulations", Feng Zhao,
  MIT Artificial Intelligence Lab, Technical Report 995, 1987.

The formula for the lexicographical ordering of a tetrahedral array was
derived by Richard Merrifield.

## Build and run the multipole program

To build and run the multipole program, type

  ./configure
  make
  multipole

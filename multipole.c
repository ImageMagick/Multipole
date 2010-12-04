/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%       M   M  U   U  L      TTTTT  IIIII  PPPP    OOO   L      EEEEE         %
%       MM MM  U   U  L        T      I    P   P  O   O  L      E             %
%       M M M  U   U  L        T      I    PPPP   O   O  L      EEE           %
%       M   M  U   U  L        T      I    P      O   O  L      E             %
%       M   M   UUU   LLLLL    T    IIIII  P       OOO   LLLLL  EEEEE         %
%                                                                             %
%                                                                             %
%          An O(N) Algorithm for Three-Dimensional N-Body Simulations         %
%                                                                             %
%                                                                             %
%                                                                             %
%                           Software Design                                   %
%                             John Cristy                                     %
%                             January 1991                                    %
%                                                                             %
%                          Algorithm Analysis                                 %
%                             Steve Lustig                                    %
%                                                                             %
%                                                                             %
%  Copyright 1999-2011 ImageMagick Studio LLC, a non-profit organization      %
%  dedicated to making software imaging solutions freely available.           %
%                                                                             %
%  You may not use this file except in compliance with the License.  You may  %
%  obtain a copy of the License at                                            %
%                                                                             %
%    http://www.imagemagick.org/script/license.php                            %
%                                                                             %
%  Unless required by applicable law or agreed to in writing, software        %
%  distributed under the License is distributed on an "AS IS" BASIS,          %
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
%  See the License for the specific language governing permissions and        %
%  limitations under the License.                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  N-body implements a sequential multipole algorithm that solves the
%  interparticle gravitational potential 1/r of N particles to within
%  round-off error and whose computational effort grows only linearly with
%  increasing N, i.e. O(N).  This makes possible simulations of a large
%  number of atoms on the order of 500,000 or more.  Typically, thousands
%  to millions of time steps are required for a realistic dynamic
%  simulation.  The linear nature of the cost of computing the potential
%  with the sequential multipole algorithm allows many more time-steps, in
%  less time, than previously possible.
%
%  The multipole-expansion algorithm computes the potential by
%  partitioning them into two parts: close-by particles (`near-field
%  potentials') and the potential of far-away particles (`far-field
%  potentials').  The near-field potential is computed using direct
%  evaluation.  The far-field potential is represented as an absolutely
%  convergent Taylor series.  Given a required precision e, the number of
%  terms p in the Taylor series can be determined in advance to obtain a
%  particular accuracy.
%
%  The goal of the multipole computation is to compute the far-field with
%  O(N) work.  This is achieved by propagating values, first upwards
%  through the computational tree, then downward in two distinct phases.
%  In the first (upward) phase the far-field interactions are computed by
%  shifting the multipole expansions of the child nodes to the center of
%  the parent node and adding the results.  The far-field interaction is
%  calculated for all cubes at all levels, beginning at the finest level
%  of refinement (the leaf nodes).  The multipole-expansion of the field
%  due to each particle is calculated relative to the center of each leaf
%  node and summed.  Proceeding upward, the expansions of the eight child
%  nodes are shifted relative to the parent's center and summed.  This
%  computation is performed at each level propagating values upward
%  through the tree.  At the end of this phase of the algorithm, each node
%  at every level has a multipole-expansion representation of the field
%  due to all particles within its boundaries.  These expansions are only
%  valid outside the region defined by the near field of each cube.
%
%  In the second (downward) phase, we form the local expansions for all
%  cubes at all levels beginning at the coarsest level (root node).  The
%  far-field multipole-expansions of all nodes in the interactive field of
%  the node under consideration are converted into local expansions
%  relative to the node's center and are summed.  Next, shift the
%  expansion of the parent cube to the center of the current cube under
%  consideration and sum with the local expansions computed from the
%  interactive field.  Note, the root node does not have a parent or
%  interactive field so it's local expansion is defined to be zero.  The
%  local expansions are propagated down the tree.  They describe the field
%  due to all particles in the system that are not contained in the
%  current node or its nearest neighbors.
%
%  Finally, for each cube i at the finest level n, we evaluate the local
%  expansion for each particle in cube i.  Each local expansion is
%  described by the coefficients of a p-term Taylor series.  Direct
%  evaluation of this series at a particle yields the potential.
%  Summating the local expansion at each leaf node generates the potential
%  or force due to all particles outside the nearest neighbor cubes. The
%  interactions of each particle in cube i are computed directly and
%  summed to generate the potential due to the nearest neightbors.
%  Finally summating far-field and near-field potential gives us the
%  desired potential.
%
%  A discussion and formal description of the multipole-expansion algorithm
%  can be found in these references:
%
%    "Computational Simulation of Multi-Body Interactions with O(N) Scaling",
%    John Cristy and David A. Pensak, Central Research and Development,
%    E.I. Du Pont de Nemours, Technical Report, June 1991.
%
%    "An O(N) Algorithm for Three-dimensional N-body Simulations", Feng Zhao,
%    MIT Artificial Intelligence Lab, Technical Report 995, 1987.
%
%  The formula for the lexicographical ordering of a tetrahedral array was
%  derived by Richard Merrifield.
%
%
*/

/*
  Include declarations.
*/
#include <multipole.h>

/*
  Global variables.
*/
static Cube
  cube;

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  B i n o m i a l                                                            %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Binomial() returns the binomial coefficient (n,k).
%
%  The format of the Binomial routine is:
%
%    Binomial(long n,long k)
%
%  A description of each parameter follows:
%
%    o n,k: Specifies an unsigned long.
%
%
*/
static unsigned long Binomial(long n,long k)
{
  register long
    i,
    j;

  unsigned long
    binomial;

  if (n < k)
    return(0);
  binomial=1;
  k=Min(k,n-k);
  j=n;
  for (i=1; i <= k; i++)
  {
    binomial=(unsigned long) ((binomial*(float) j/i)+0.5);
    j--;
  }
  return(binomial);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  C l a s s i f i c a t i o n                                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Classification() creates a linked list of particles in the leaf nodes of the
%  particle cube.
%
%  The format of the Classification routine is:
%
%    Classification(Particle *particles,unsigned long number_particles)
%
%  A description of each parameter follows:
%
%    o particles: Specifies a pointer to an array of Particles structures.
%
%    o number_particles: The number of particules in the Particles array.
%
*/
static void Classification(Particle *particles,unsigned long number_particles)
{
  register long
    i,
    level;

  register Node
    *node;

  register Particle
    *p;

  unsigned long
    id;

  p=particles;
  for (i=0; i < (long) number_particles; i++)
  {
    /*
      Descend the tree to a particular leaf node.
    */
    node=cube.root;
    for (level=1; level < (long) cube.depth; level++)
    {
      id=(p->x > node->mid_x) | (p->y > node->mid_y) << 1 |
        (p->z > node->mid_z) << 2;
      node=node->child[id];
    }
    p->next=node->particle;
    node->particle=p;
    p++;
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  C o m p u t e P h i                                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ComputePhi() forms the multipole expansions of the potential field due to
%  particle in each node.
%
%  The format of the ComputePhi routine is:
%
%    ComputePhi(register Node *node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void ComputePhi(register Node *node)
{
  double
    x_distance,
    y_distance,
    z_distance;

  long
    i,
    j,
    k;

  register double
    sum;

  register float
    *ijk,
    *phi;

  register long
    a,
    b,
    g;

  unsigned long
    id;

  /*
    Recursively descend the particle cube.
  */
  if (node->level < (cube.depth-1))
    for (id=0; id < MaxChildren; id++)
      ComputePhi(node->child[id]);
  if (node->level == (cube.depth-1))
    if (node->particle != (Particle *) NULL)
      {
        register Particle
          *p;

        /*
          Form multipole expansions of potential field due to particles
          in each cube about the cube center at the leaf nodes.
        */
        for (p=node->particle; p != (Particle *) NULL; p=p->next)
        {
          x_distance=p->x-node->mid_x;
          y_distance=p->y-node->mid_y;
          z_distance=p->z-node->mid_z;
          for (i=0; i <= (long) cube.precision; i++)
          {
            cube.x_power[i]=pow(x_distance,(double) i);
            cube.y_power[i]=pow(y_distance,(double) i);
            cube.z_power[i]=pow(z_distance,(double) i);
          }
          phi=node->phi;
          ijk=cube.ijk_factorial;
          for (i=0; i <= (long) cube.precision; i++)
            for (j=0; j <= (long) (cube.precision-i); j++)
              for (k=0; k <= (long) (cube.precision-i-j); k++)
              {
                sum=cube.x_power[i]*cube.y_power[j]*cube.z_power[k]*
                  cube.one_power[i+j+k]*(*ijk++);
                *phi+=sum;
                phi++;
              }
        }
      }
  if (node->level > 1)
    {
      /*
        Form multipole expansions about the centers of each cube at all
        internal nodes, each expansion representing the potential field
        due to all particles contained in the cube.  Shift the center
        of the cube's expansion to the parent cube center using the
        `Translation of a Multipole Expansion' lemma.
      */
      MultipoleExpansion(node,node->parent,cube.phi_tilde);
      phi=node->parent->phi;
      for (i=0; i <= (long) cube.precision; i++)
        for (j=0; j <= (long) (cube.precision-i); j++)
          for (k=0; k <= (long) (cube.precision-i-j); k++)
          {
            sum=0.0;
            for (a=0; a <= i; a++)
              for (b=0; b <= j; b++)
                for (g=0; g <= k; g++)
                  sum+=TetrahedralArray(node->phi,a,b,g)*
                    TetrahedralArray(cube.phi_tilde,i-a,j-b,k-g);
            *phi+=sum;
            phi++;
          }
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  C o m p u t e P s i                                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ComputePsi() forms a local expansion which describes the potential field due
%  to all particles in the system that are not contained in the current cube or
%  its near-field.
%
%  The format of the ComputePsi routine is:
%
%    ComputePsi(register Node *node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void ComputePsi(register Node *node)
{
  double
    x_distance,
    y_distance,
    z_distance;

  long
    i,
    j,
    k,
    n;

  Node
    *parent;

  register long
    a,
    b,
    g;

  register double
    sum;

  register float
    *ijk,
    *phi,
    *psi;

  register Node
    **interactive_neighbor;

  unsigned long
    id;

  if (node->level == (cube.depth-1))
    {
      /*
        Form a local expansion by using the `Conversion of Multipole
        into a Local Expansion' lemma to convert the multipole expansion
        of each cube in the interaction-field of the current leaf cube
        to a local expansion about the center of the leaf.
      */
      interactive_neighbor=node->interactive_field;
      while (*interactive_neighbor != (Node *) NULL)
      {
        /*
          Shift multipole expansion of the interactive cube to the
          center of the current leaf cube.
        */
        LocalExpansion(*interactive_neighbor,node,cube.psi_tilde);
        psi=node->psi;
        phi=(*interactive_neighbor)->phi;
        ijk=cube.ijk_factorial;
        for (i=0; i <= (long) cube.precision; i++)
          for (j=0; j <= (long) (cube.precision-i); j++)
            for (k=0; k <= (long) (cube.precision-i-j); k++)
            {
              sum=0.0;
              n=i+j+k;
              for (a=0; a <= (long) (cube.precision-n); a++)
                for (b=0; b <= (long) (cube.precision-n-a); b++)
                  for (g=0; g <= (long) (cube.precision-n-a-b); g++)
                    sum+=TetrahedralArray(phi,a,b,g)*
                      TetrahedralArray(cube.psi_tilde,i+a,j+b,k+g);
              *psi+=sum*(*ijk++);
              psi++;
            }
        interactive_neighbor++;
      }
      /*
        Shift local expansion of the parent cube to the center of the
        current leaf cube using the `Translation of a Local Expansion'
        lemma.
      */
      parent=node->parent;
      x_distance=node->mid_x-parent->mid_x;
      y_distance=node->mid_y-parent->mid_y;
      z_distance=node->mid_z-parent->mid_z;
      for (i=0; i <= (long) cube.precision; i++)
      {
        cube.x_power[i]=pow(x_distance,(double) i);
        cube.y_power[i]=pow(y_distance,(double) i);
        cube.z_power[i]=pow(z_distance,(double) i);
      }
      psi=node->psi;
      ijk=cube.ijk_factorial;
      for (i=0; i <= (long) cube.precision; i++)
        for (j=0; j <= (long) (cube.precision-i); j++)
          for (k=0; k <= (long) (cube.precision-i-j); k++)
          {
            sum=0.0;
            n=i+j+k;
            for (a=0; a <= (long) (cube.precision-n); a++)
              for (b=0; b <= (long) (cube.precision-n-a); b++)
                for (g=0; g <= (long) (cube.precision-n-a-b); g++)
                  sum+=TetrahedralArray(parent->psi,i+a,j+b,k+g)*
                    TetrahedralArray(cube.ijk_factorial,a,b,g)*
                    cube.x_power[a]*cube.y_power[b]*cube.z_power[g];
            *psi+=sum*(*ijk++);
            psi++;
          }
      EvaluatePotential(node);
    }
  else
    if (node->level > 1)
      {
        /*
          Form a local expansion by using the `Conversion of Multipole
          into a Local Expansion' lemma to convert the multipole
          expansion of each cube in the interaction-field of the
          current cube to a local expansion about the center of the
          cube.
        */
        interactive_neighbor=node->interactive_field;
        while (*interactive_neighbor != (Node *) NULL)
        {
          /*
            Shift multipole expansion of the interactive cube to the
            center of the current cube.
          */
          LocalExpansion(*interactive_neighbor,node,cube.psi_tilde);
          psi=node->psi;
          phi=(*interactive_neighbor)->phi;
          for (i=0; i <= (long) cube.precision; i++)
            for (j=0; j <= (long) (cube.precision-i); j++)
              for (k=0; k <= (long) (cube.precision-i-j); k++)
              {
                sum=0.0;
                n=i+j+k;
                for (a=0; a <= (long) (cube.precision-n); a++)
                  for (b=0; b <= (long) (cube.precision-n-a); b++)
                    for (g=0; g <= (long) (cube.precision-n-a-b); g++)
                      sum+=TetrahedralArray(phi,a,b,g)*
                        TetrahedralArray(cube.psi_tilde,i+a,j+b,k+g);
                *psi+=sum;
                psi++;
              }
          interactive_neighbor++;
        }
        /*
          Shift local expansion of the parent cube to the center of the
          current node using the `Translation of a Local Expansion'.
        */
        parent=node->parent;
        x_distance=node->mid_x-parent->mid_x;
        y_distance=node->mid_y-parent->mid_y;
        z_distance=node->mid_z-parent->mid_z;
        for (i=0; i <= (long) cube.precision; i++)
        {
          cube.x_power[i]=pow(x_distance,(double) i);
          cube.y_power[i]=pow(y_distance,(double) i);
          cube.z_power[i]=pow(z_distance,(double) i);
        }
        psi=node->psi;
        for (i=0; i <= (long) cube.precision; i++)
          for (j=0; j <= (long) (cube.precision-i); j++)
            for (k=0; k <= (long) (cube.precision-i-j); k++)
            {
              sum=0.0;
              n=i+j+k;
              for (a=0; a <= (long) (cube.precision-n); a++)
                for (b=0; b <= (long) (cube.precision-n-a); b++)
                  for (g=0; g <= (long) (cube.precision-n-a-b); g++)
                    sum+=TetrahedralArray(parent->psi,i+a,j+b,k+g)*
                      TetrahedralArray(cube.ijk_factorial,a,b,g)*
                      cube.x_power[a]*cube.y_power[b]*cube.z_power[g];
              *psi+=sum;
              psi++;
            }
      }
  /*
    Recursively descend the particle cube.
  */
  if (node->level < (cube.depth-1))
    for (id=0; id < MaxChildren; id++)
      ComputePsi(node->child[id]);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  C r e a t e C u b e                                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CreateCube() creates all levels of the particle cube.
%
%  The format of the CreateCube routine is:
%
%    CreateCube(register Node *node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void CreateCube(register Node *node)
{
  register double
    bisect;

  register unsigned long
    id,
    level;

  /*
    Create a child cube and recursively descend to the next level.
  */
  level=node->level+1;
  bisect=cube.edge_length[level]*0.5;
  for (id=0; id < MaxChildren; id++)
  {
    node->child[id]=InitializeNode(id,level,node,
      node->mid_x+(id & 1 ? bisect : -bisect),
      node->mid_y+(id & 2 ? bisect : -bisect),
      node->mid_z+(id & 4 ? bisect : -bisect));
    if (node->child[id] == (Node *) NULL)
      Error("unable to allocate memory",(char *) NULL);
    if (level < (cube.depth-1))
      CreateCube(node->child[id]);
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  D e f i n e I n t e r a c t i v e F i e l d                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DefineInteractiveField() assigns the near-field of a particular cube.  The
%  interaction-field of a node is the set of nodes which are children of the
%  nodes in the near-field of it parent and which are not in the near-field of
%  itself.  If all eight children of a parent node is in the interaction-field,
%  they are replaced by their common parent node (super_node).
%
%  The format of the DefineInteractiveField routine is:
%
%    DefineInteractiveField(node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void DefineInteractiveField(register Node *node)
{
  register Node
    **interactive_neighbor,
    **near_neighbor;

  register unsigned long
    id,
    super_node;

  /*
    Recursively descend the particle cube.
  */
  if (node->level < (cube.depth-1))
    for (id=0; id < MaxChildren; id++)
      DefineInteractiveField(node->child[id]);
  if (node->level > 0)
    {
      /*
        Nodes in the parents near-field are candidates for the interactive
        field.
      */
      interactive_neighbor=node->interactive_field;
      near_neighbor=node->parent->near_field;
      while (*near_neighbor != (Node *) NULL)
      {
        super_node=True;
        for (id=0; id < MaxChildren; id++)
          if (InNearField((*near_neighbor)->child[id],node))
            super_node=False;
          else
            {
              *interactive_neighbor=(*near_neighbor)->child[id];
              interactive_neighbor++;
            }
        super_node=False;  /* for our benchmark, turn supernode off */
        if (super_node)
          {
            /*
              Supernode:  All 8 children represented by this node.
            */
            interactive_neighbor-=MaxChildren;
            *interactive_neighbor=(*near_neighbor);
            interactive_neighbor++;
          }
        near_neighbor++;
      }
      *interactive_neighbor=(Node *) NULL;
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  D e f i n e N e a r F i e l d                                              %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DefineNearField() assigns the near-field of a particular cube.  The
%  near-field of a node is the set of nodes at the same level of refinement
%  which shares a boundary point with it.
%
%  The format of the DefineNearField routine is:
%
%    DefineNearField(node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void DefineNearField(register Node *node)
{
  register unsigned long
    id;

  /*
    Recursively descend the particle cube.
  */
  if (node->level < (cube.depth-1))
    for (id=0; id < MaxChildren; id++)
      DefineNearField(node->child[id]);
  if (node->level > 0)
    {
      /*
        Find the near neighbors of this particular node.
      */
      cube.near_neighbor=node->near_field;
      FindNeighbors(cube.root,node);
      *cube.near_neighbor=(Node *) NULL;
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   E r r o r                                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Error() displays an error message and then terminates the program.
%
%  The format of the Error routine is:
%
%      Error(char *message,char *qualifier)
%
%  A description of each parameter follows:
%
%    o message: Specifies the message to display before terminating the
%      program.
%
%    o qualifier: Specifies any qualifier to the message.
%
*/
static void Error(char *message,char *qualifier)
{
  (void) fprintf(stderr,"%s: %s",application_name,message);
  if (qualifier != (char *) NULL)
    (void) fprintf(stderr," (%s)",qualifier);
  (void) fprintf(stderr,".\n");
  exit(1);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  E v a l u a t e P o t e n t i a l                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  EvaluatePotential() evaluates the local expansions at the particle position
%  to determine the far field potential.  The direct near-field potential is
%  summed to give the total potential.
%
%  The format of the EvaluatePotential routine is:
%
%    EvaluatePotential(register Node *node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static void EvaluatePotential(register Node *node)
{
  double
    distance_squared,
    far_potential,
    near_potential,
    x_distance,
    y_distance,
    z_distance;

  register float
    *psi;

  register long
    i,
    j,
    k;

  register Node
    **near_neighbor;

  register Particle
    *p,
    *q;

  far_potential=0.0;
  near_potential=0.0;
  for (p=node->particle; p != (Particle *) NULL; p=p->next)
  {
    /*
      Evaluate local expansions at particle positions.
    */
    x_distance=p->x-node->mid_x;
    y_distance=p->y-node->mid_y;
    z_distance=p->z-node->mid_z;
    for (i=0; i <= (long) cube.precision; i++)
    {
      cube.x_power[i]=pow(x_distance,(double) i);
      cube.y_power[i]=pow(y_distance,(double) i);
      cube.z_power[i]=pow(z_distance,(double) i);
    }
    psi=node->psi;
    for (i=0; i <= (long) cube.precision; i++)
      for (j=0; j <= (long) (cube.precision-i); j++)
        for (k=0; k <= (long) (cube.precision-i-j); k++)
        {
          far_potential+=(*psi)*cube.x_power[i]*cube.y_power[j]*cube.z_power[k];
          psi++;
        }
    /*
      Compute potential due to particles in current cube directly.  Use
      Newton's third law to reduce the number of pairwise interactions.
    */
    for (q=p->next; q != (Particle *) NULL; q=q->next)
    {
      x_distance=p->x-q->x;
      y_distance=p->y-q->y;
      z_distance=p->z-q->z;
      distance_squared=x_distance*x_distance+y_distance*y_distance+
        z_distance*z_distance;
      near_potential+=1.0/sqrt(distance_squared);
    }
    /*
      Compute potential due to particles in near-field directly.
    */
    near_neighbor=node->near_field;
    while (*near_neighbor != (Node *) NULL)
    {
      for (q=(*near_neighbor)->particle; q != (Particle *) NULL; q=q->next)
      {
        if (p > q)
          continue;  /* canonical form */
        x_distance=p->x-q->x;
        y_distance=p->y-q->y;
        z_distance=p->z-q->z;
        distance_squared=x_distance*x_distance+y_distance*y_distance+
          z_distance*z_distance;
        near_potential+=1.0/sqrt(distance_squared);
      }
      near_neighbor++;
    }
  }
  cube.near_potential+=2.0*GravitationalConstant*near_potential;
  cube.far_potential+=GravitationalConstant*far_potential;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  F i n d N e i g h b o r s                                                  %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FindNeighbors() determines which nodes in the particle cube are in the near-
%  field of the reference node.  The near-field of a cube is defined as
%  consisiting of those cubes that are on the same level as the cube, and have
%  a distance to the center of the cube less than sqrt(3) of the cube edge
%  length.
%
%  The format of the FindNeighbors routine is:
%
%    FindNeighbors(register Node *node,register Node *reference_node)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
%    o reference_node: Specifies a pointer to a Node structure.
%
*/
static void FindNeighbors(register Node *node,register Node *reference_node)
{
  register double
    distance_squared;

  register double
    x_distance,
    y_distance,
    z_distance;

  unsigned long
    id;

  /*
    Recursively descend the particle cube.
  */
  if (node->level < reference_node->level)
    for (id=0; id < MaxChildren; id++)
      FindNeighbors(node->child[id],reference_node);
  if (node->level == reference_node->level)
    if (node != reference_node)
      {
        /*
          Determine if this node is within the specified distance of the
          reference node.
        */
        x_distance=reference_node->mid_x-node->mid_x;
        y_distance=reference_node->mid_y-node->mid_y;
        z_distance=reference_node->mid_z-node->mid_z;
        distance_squared=x_distance*x_distance+y_distance*y_distance+
          z_distance*z_distance;
        if (distance_squared <= cube.diagonal_length[node->level])
          {
            *cube.near_neighbor=node;
            cube.near_neighbor++;
          }
      }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  I n i t i a l i z e N o d e                                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  InitializeNode() allocates memory for a new node in the color cube treee and
%  presets all fields to zero.
%
%  The format of the InitializeNode routine is:
%
%      node=InitializeNode(unsigned long id,unsigned long level,Node *parent,
%        double mid_x,double mid_y,double mid_z)
%
%  A description of each parameter follows.
%
%    o node: The InitializeNode routine returns this integer address.
%
%    o id: Specifies the child number of the node.
%
%    o level: Specifies the level in the particle cube the node resides.
%
%    o parent: Specifies the parent of the initialized node.
%
%    o mid_x: Specifies the mid point of the x extent for this node.
%
%    o mid_y: Specifies the mid point of the y extent for this node.
%
%    o mid_z: Specifies the mid point of the z extent for this node.
%
*/
static Node *InitializeNode(unsigned long id,unsigned long level,Node *parent,
  double mid_x,double mid_y,double mid_z)
{
  register float
    *phi,
    *psi;

  register long
    i;

  register Node
    *node;

  if (cube.free_nodes == 0)
    {
      Nodes
        *nodes;

      /*
        Allocate a queue of nodes.
      */
      nodes=(Nodes *) malloc(sizeof(Nodes));
      if (nodes == (Nodes *) NULL)
        return((Node *) NULL);
      nodes->phi=(float *)
        malloc(NodesInAQueue*cube.coefficients*sizeof(float));
      if (nodes->phi == (float *) NULL)
        return((Node *) NULL);
      nodes->psi=(float *)
        malloc(NodesInAQueue*cube.coefficients*sizeof(float));
      if (nodes->phi == (float *) NULL)
        return((Node *) NULL);
      nodes->next=cube.node_queue;
      cube.node_queue=nodes;
      cube.next_node=nodes->nodes;
      cube.free_nodes=NodesInAQueue;
      cube.phi=nodes->phi;
      cube.psi=nodes->psi;
    }
  /*
    Initialize node fields.
  */
  node=cube.next_node;
  node->parent=parent;
  for (i=0; i < MaxChildren; i++)
    node->child[i]=(Node *) NULL;
  node->near_field[0]=(Node *) NULL;
  node->interactive_field[0]=(Node *) NULL;
  node->id=id;
  node->level=level;
  node->mid_x=mid_x;
  node->mid_y=mid_y;
  node->mid_z=mid_z;
  node->phi=cube.phi;
  node->psi=cube.psi;
  /*
    Set multipole and local expansions to zero.
  */
  phi=node->phi;
  psi=node->psi;
  for (i=0; i < (long) cube.coefficients; i++)
  {
    *phi++=0.0;
    *psi++=0.0;
  }
  node->particle=(Particle *) NULL;
  /*
    Prepare for next node.
  */
  cube.nodes++;
  cube.free_nodes--;
  cube.next_node++;
  cube.phi+=cube.coefficients;
  cube.psi+=cube.coefficients;
  return(node);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  I n N e a r F i e l d                                                      %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  InNearField() returns True if the specified node is a member of the
%  near-field.
%
%  The format of the InNearField routine is:
%
%    InNearField(register Node *reference_node,register Node *node)
%
%  A description of each parameter follows:
%
%    o reference_node: Specifies a pointer to a Node structure.
%
%    o node: Specifies a pointer to a Node structure.
%
*/
static unsigned long InNearField(register Node *reference_node,
  register Node *node)
{
  register Node
    **near_neighbor;

  if (reference_node == node)
    return(True);
  /*
    Search the near-field for the specified node.
  */
  near_neighbor=node->near_field;
  while (*near_neighbor != (Node *) NULL)
  {
    if (*near_neighbor == reference_node)
      return(True);
    near_neighbor++;
  }
  return(False);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  L o c a l E x p a n s i o n                                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LocalExpansion() forms a local expansion which describes the potential field
%  due to all particles in the system that are not contained in the current
%  cube or its near-field.
%
%  The format of the LocalExpansion routine is:
%
%    LocalExpansion(Node *interactive_neighbor,Node *node,float *psi_tilde)
%
%  A description of each parameter follows:
%
%    o interactive_neighbor: Specifies a pointer to a Node structure.
%
%    o node: Specifies a pointer to a Node structure.
%
%    o psi_tilde: Specifies a pointer to an array of floats representing the
%      local expansion for this node.
%
*/
static void LocalExpansion(Node *interactive_neighbor,Node *node,
 float *psi_tilde)
{
  double
    distance,
    tau_x,
    tau_y,
    tau_z,
    x_distance,
    y_distance,
    z_distance;

  long
    i,
    j,
    k,
    n;

  register double
    sum;

  register float
    *ijk;

  register long
    a,
    b,
    g;

  static float
    r_zero[HighestDegreeHarmonic+1];

  /*
    Initialize Psi' as specified by Theorem 3.2.2.
  */
  x_distance=interactive_neighbor->mid_x-node->mid_x;
  y_distance=interactive_neighbor->mid_y-node->mid_y;
  z_distance=interactive_neighbor->mid_z-node->mid_z;
  distance=sqrt(x_distance*x_distance+y_distance*y_distance+
    z_distance*z_distance);
  tau_x=2.0*(x_distance/distance);
  tau_y=2.0*(y_distance/distance);
  tau_z=2.0*(z_distance/distance);
  for (i=0; i <= (long) cube.precision; i++)
  {
    cube.x_power[i]=pow(tau_x,(double) i);
    cube.y_power[i]=pow(tau_y,(double) i);
    cube.z_power[i]=pow(tau_z,(double) i);
    r_zero[i]=1.0/pow(distance,(double) i+1);
  }
  ijk=cube.ijk_factorial;
  for (i=0; i <= (long) cube.precision; i++)
    for (j=0; j <= (long) (cube.precision-i); j++)
      for (k=0; k <= (long) (cube.precision-i-j); k++)
      {
        sum=0.0;
        for (a=0; a <= (i/2); a++)
          for (b=0; b <= (j/2); b++)
            for (g=0; g <= (k/2); g++)
              sum+=TetrahedralArray(cube.ijk_binomial,i-a,j-b,k-g)*
                cube.binomial[i-a][a]*cube.binomial[j-b][b]*
                cube.binomial[k-g][g]*cube.x_power[i-2*a]*cube.y_power[j-2*b]*
                cube.z_power[k-2*g];
        n=i+j+k;
        sum*=cube.one_power[n]*r_zero[n]/(*ijk++);
        *psi_tilde++=sum;
      }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M u l t i p o l e E x p a n s i o n                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MultipoleExpansion() forms the multipole expansions of the potential field
%  due to particle in each node.
%
%  The format of the MultipoleExpansion routine is:
%
%    MultipoleExpansion(Node *node,Node *parent,register float *phi_tilde)
%
%  A description of each parameter follows:
%
%    o node: Specifies a pointer to a Node structure.
%
%    o parent: Specifies a pointer to a Node structure.
%
%    o phi_tilde: Specifies a pointer to an array of floats representing the
%      multipole expansion for this node.
%
*/
static void MultipoleExpansion(Node *node,Node *parent,
  register float *phi_tilde)
{
  double
    x_distance,
    y_distance,
    z_distance;

  register float
    *ijk;

  register long
    i,
    j,
    k;

  /*
    Initialize Phi' as specified by Theorem 3.2.1.
  */
  x_distance=node->mid_x-parent->mid_x;
  y_distance=node->mid_y-parent->mid_y;
  z_distance=node->mid_z-parent->mid_z;
  for (i=0; i <= (long) cube.precision; i++)
  {
    cube.x_power[i]=pow(-x_distance,(double) i);
    cube.y_power[i]=pow(-y_distance,(double) i);
    cube.z_power[i]=pow(-z_distance,(double) i);
  }
  ijk=cube.ijk_factorial;
  for (i=0; i <= (long) cube.precision; i++)
    for (j=0; j <= (long) (cube.precision-i); j++)
      for (k=0; k <= (long) (cube.precision-i-j); k++)
        *phi_tilde++=cube.x_power[i]*cube.y_power[j]*cube.z_power[k]*(*ijk++);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M u l t i p o l e P o t e n t i a l                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MultipolePotential() calculates the interaction potentials of a number of
%  particles in O(N).
%
%  The format of the MultipolePotential routine is:
%
%    potential=MultipolePotential(Particle *particles,
%      unsigned long number_particles,unsigned long precision,
%      unsigned long tree_depth,double minimum_extent,double maximum_extent)
%
%  A description of each parameter follows:
%
%    o particles: Specifies a pointer to an array of Particles structures.
%
%    o number_particles: Specifies the number of particles in the
%      particles array.
%
%    o precision: Specifies the number of terms to use in the multipole
%      and local expansions.
%
%    o tree_depth: Normally, this integer value is zero or one.  A zero
%      tells MultipolePotential to choose a optimal tree depth of
%      Log8(number_particles).
%
%    o minimum_extent: Specifies the minumum of the cube extent.
%
%    o maximum_extent: Specifies the maximum of the cube extent.
%
*/
static double MultipolePotential(Particle *particles,
  unsigned long number_particles,unsigned long precision,
  unsigned long tree_depth,double minimum_extent,double maximum_extent)
{
  double
    factorial[HighestDegreeHarmonic+1];

  long
    count,
    start_time;

  Nodes
    *nodes;

  register double
    sum;

  register float
    *ijk,
    *ijk_binomial;

  register long
    i,
    j,
    k,
    n;

  unsigned long
    level;

  (void) fprintf(stderr,"Particles: %lu\n",number_particles);
  cube.precision=precision;
  if (cube.precision > HighestDegreeHarmonic)
    cube.precision=HighestDegreeHarmonic;
  (void) fprintf(stderr,"Precision: %lu\n",cube.precision);
  cube.coefficients=Binomial((int) cube.precision+3,3);
  (void) fprintf(stderr,"Coefficients: %lu\n",cube.coefficients);
  /*
    Initialize ijk factorial table.
  */
  factorial[0]=1.0;
  sum=1.0;
  for (i=1; i <= (long) cube.precision; i++)
  {
    sum*=(double) i;
    factorial[i]=sum;
  }
  cube.ijk_factorial=(float *) malloc(cube.coefficients*sizeof(float));
  if (cube.ijk_factorial == (float *) NULL)
    Error("unable to allocate memory",(char *) NULL);
  ijk=cube.ijk_factorial;
  for (i=0; i <= (long) cube.precision; i++)
    for (j=0; j <= (long) (cube.precision-i); j++)
      for (k=0; k <= (long) (cube.precision-i-j); k++)
        *ijk++=1.0/(factorial[i]*factorial[j]*factorial[k]);
  /*
    Initialize ijk binomial table: [(-1/2)*(-1/2-1)*(-1/2-2)*...].
  */
  factorial[0]=1.0;
  sum=1.0;
  for (i=1; i <= (long) cube.precision; i++)
  {
    sum*=(-0.5-(double) i+1.0);
    factorial[i]=sum;
  }
  cube.ijk_binomial=(float *) malloc(cube.coefficients*sizeof(float));
  if (cube.ijk_binomial == (float *) NULL)
    Error("unable to allocate memory",(char *) NULL);
  ijk_binomial=cube.ijk_binomial;
  ijk=cube.ijk_factorial;
  for (i=0; i <= (long) cube.precision; i++)
    for (j=0; j <= (long) (cube.precision-i); j++)
      for (k=0; k <= (long) (cube.precision-i-j); k++)
        *ijk_binomial++=factorial[i+j+k]*(*ijk++);
  /*
    Initialize binomial coefficient table.
  */
  for (n=0; n <= (long) cube.precision; n++)
    for (k=0; k <= (long) cube.precision; k++)
      cube.binomial[n][k]=(float) Binomial(n,k);
  /*
    Initialize power of one lookup table.
  */
  for (i=0; i <= (long) cube.precision; i++)
    cube.one_power[i]=pow((double) -1.0,(double) i);
  /*
    Initialize Tetrahedral cubic array lookup table.
  */
  for (i=0; i <= (long) (cube.precision+2); i++)
    cube.tetra_three[i]=cube.coefficients-Binomial((int) cube.precision-i+2,3);
  for (i=0; i <= (long) (cube.precision+2); i++)
    cube.tetra_two[i]=Binomial((int) cube.precision-i+2,2);
  /*
    Allocate Phi' & Psi'.
  */
  cube.phi_tilde=(float *) malloc(cube.coefficients*sizeof(float));
  cube.psi_tilde=(float *) malloc(cube.coefficients*sizeof(float));
  if ((cube.phi_tilde == (float *) NULL) ||
      (cube.psi_tilde == (float *) NULL))
    Error("unable to allocate memory",(char *) NULL);
  cube.near_potential=0.0;
  cube.far_potential=0.0;
  /*
    Choose a level of refinement:  depth is log8(number_particles)-1.
  */
  if (tree_depth > 1)
    cube.depth=Min(tree_depth,8);
  else
    {
      cube.depth=0;
      count=(long) (number_particles >> 3);
      do
      {
        count>>=3;
        cube.depth++;
      }
      while (count > 0);
      if (cube.depth < 1)
        cube.depth=1;
    }
  (void) fprintf(stderr,"Cube depth: %lu\n",cube.depth);
  /*
    Initialize edge length table.
  */
  (void) fprintf(stderr,"Cube lower left front vertex: %.8g %.8g %.8g\n",
    minimum_extent,minimum_extent,minimum_extent);
  (void) fprintf(stderr,"Cube upper right behind vertex: %.8g %.8g %.8g\n",
    maximum_extent,maximum_extent,maximum_extent);
  cube.edge_length[0]=maximum_extent-minimum_extent;
  cube.diagonal_length[0]=2.0*CircumscribedSphereRadius*cube.edge_length[0]*
    cube.edge_length[0];
  for (level=1; level < cube.depth; level++)
  {
    cube.edge_length[level]=cube.edge_length[level-1]*0.5;
    cube.diagonal_length[level]=2.0*CircumscribedSphereRadius*
      cube.edge_length[level]*cube.edge_length[level];
  }
  /*
    Allocate particle description tree.
  */
  cube.node_queue=(Nodes *) NULL;
  cube.nodes=0;
  cube.free_nodes=0;
  cube.root=InitializeNode(0,0,(Node *) NULL,
    minimum_extent+cube.edge_length[0]*0.5,
    minimum_extent+cube.edge_length[0]*0.5,
    minimum_extent+cube.edge_length[0]*0.5);
  if (cube.root == (Node *) NULL)
    Error("unable to allocate memory",(char *) NULL);
  cube.root->parent=cube.root;
  /*
    Multibody steps:
  */
  (void) fprintf(stderr,"CreateCube: ");
  start_time=ProcessTime();
  CreateCube(cube.root);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"Classification: ");
  start_time=ProcessTime();
  Classification(particles,number_particles);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"DefineNearField: ");
  start_time=ProcessTime();
  DefineNearField(cube.root);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"DefineInteractiveField: ");
  start_time=ProcessTime();
  DefineInteractiveField(cube.root);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"ComputePhi: ");
  start_time=ProcessTime();
  ComputePhi(cube.root);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"ComputePsi: ");
  start_time=ProcessTime();
  ComputePsi(cube.root);
  (void) fprintf(stderr,"%lds\n",ProcessTime()-start_time);
  (void) fprintf(stderr,"Near: %.8g\nFar: %.8g\n",cube.near_potential,
    cube.far_potential);
  /*
    Release N-body cube tree storage.
  */
  do
  {
    nodes=cube.node_queue->next;
    (void) free((char *) cube.node_queue->phi);
    (void) free((char *) cube.node_queue->psi);
    (void) free((char *) cube.node_queue);
    cube.node_queue=nodes;
  }
  while (cube.node_queue != (Nodes *) NULL);
  (void) free((char *) cube.ijk_factorial);
  (void) free((char *) cube.ijk_binomial);
  (void) free((char *) cube.phi_tilde);
  (void) free((char *) cube.psi_tilde);
  return(cube.near_potential+cube.far_potential);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  N a i v e P o t e n t i a l                                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  NaivePotential() calculates the interaction potentials of a number of
%  particles using the naive direct method.
%
%  The format of the NaivePotential routine is:
%
%    potential=NaivePotential(Particle *particles,
%      long number_particles)
%
%  A description of each parameter follows:
%
%    o particles: Specifies a pointer to an array of Particles structures.
%
%    o number_particles: Specifies the number of particles in the
%      particles array.
%
*/
static double NaivePotential(Particle *particles,unsigned long number_particles)
{
  double
    distance_squared,
    total_potential,
    x_distance,
    y_distance,
    z_distance;

  register int
    i,
    j;

  register Particle
    *p,
    *q;

  /*
    Calculate gravitational total potential.  Use Newton's third law to
    reduce the number of pairwise interactions.
  */
  total_potential=0.0;
  p=particles;
  for (i=0; i < (long) number_particles; i++)
  {
    q=p;
    for (j=(i+1); j < (long) number_particles; j++)
    {
      q++;
      x_distance=p->x-q->x;
      y_distance=p->y-q->y;
      z_distance=p->z-q->z;
      distance_squared=x_distance*x_distance+y_distance*y_distance+
        z_distance*z_distance;
      total_potential+=1.0/sqrt(distance_squared);
    }
    p++;
  }
  return(2.0*GravitationalConstant*total_potential);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   P r o c e s s T i m e                                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ProcessTime() returns the number of processor seconds our program has
%  consumed.
%
%  The format of the ProcessTime routine is:
%
%      seconds=ProcessTime(void)
%
%  A description of each parameter follows:
%
%    o seconds: return the number of processor seconds our program has consumed.
%
*/
static long ProcessTime(void)
{
#if defined(mac)
   return((long) time(0));
#elif defined(WIN32)
   return(GetTickCount()/1000);
#else
#include <sys/times.h>
#include <limits.h>

#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

  struct tms
    usage;

  (void) times(&usage);
  return((long) (usage.tms_utime/CLK_TCK));
#endif
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   M a i n                                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
*/
int main(argc,argv)
int
  argc;

char
  **argv;
{
#define MinRange  0.0
#define MaxRange  255.9999999

  double
    multipole_potential,
    maximum_extent,
    minimum_extent,
    naive_potential;

  long
    start_time;

  Particle
    *particles;

  register long
    i;

  unsigned long
    number_particles,
    precision;

  /*
    Allocate particle array.
  */
  application_name=argv[0];
  number_particles=4095;
  if (argc > 1)
    number_particles=(unsigned long) atol(argv[1]);
  precision=3;
  if (argc > 2)
    precision=(unsigned long) atol(argv[2]);
  particles=(Particle *) malloc(number_particles*sizeof(Particle));
  if (particles == (Particle *) NULL)
    Error("unable to allocate memory",(char *) NULL);
  /*
    Read particles positions.
  */
  srand(30428877);  /* must be fixed to obtain proper benchmarked results */
  minimum_extent=MaxRange;
  maximum_extent=MinRange;
  for (i=0; i < (long) number_particles; i++)
  {
    /*
      Define vertexes of the particle cube.
    */
    particles[i].x=(rand()*((MaxRange-MinRange)+MinRange))/RAND_MAX;
    particles[i].y=(rand()*((MaxRange-MinRange)+MinRange))/RAND_MAX;
    particles[i].z=(rand()*((MaxRange-MinRange)+MinRange))/RAND_MAX;
    if (particles[i].x < minimum_extent)
      minimum_extent=particles[i].x;
    if (particles[i].y < minimum_extent)
      minimum_extent=particles[i].y;
    if (particles[i].z < minimum_extent)
      minimum_extent=particles[i].z;
    if (particles[i].x > maximum_extent)
      maximum_extent=particles[i].x;
    if (particles[i].y > maximum_extent)
      maximum_extent=particles[i].y;
    if (particles[i].z > maximum_extent)
      maximum_extent=particles[i].z;
  }
  minimum_extent=floor(minimum_extent);
  maximum_extent=ceil(maximum_extent);
  start_time=ProcessTime();
  multipole_potential=MultipolePotential(particles,number_particles,precision,
    0,minimum_extent,maximum_extent);
  (void) fprintf(stderr,"The multipole-expansion potential is %.8g (%lds).\n",
    multipole_potential,ProcessTime()-start_time);
  start_time=ProcessTime();
  naive_potential=NaivePotential(particles,number_particles);
  (void) fprintf(stderr,"The naive potential is %.8g (%lds).\n",
    naive_potential,ProcessTime()-start_time);
  (void) fprintf(stderr,"The expansion error: %.8g\n\n",
    AbsoluteValue(multipole_potential-naive_potential));
  (void) free((char *) particles);
  return(False);
}

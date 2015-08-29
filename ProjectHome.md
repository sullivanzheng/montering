## A monte carlo simulation project for DNA modeling ##

**General Description**

Could DNA local structure influence its global topology and the topological selectivity of the recombination product of serine recombinase? This remains an unknown topic.

We investigate this problem by adding local rigid body model inferred from crystal structure  to worm-like chain modeled DNA rings and, using Metropolis Monte Carlo simulation and various sampling methods, including Gibbs sampling, umbrella sampling or other biased sampling methods, the rare events of DNA _**res**_ site synapsis with distinct topologies are observed and sampled. Their conditional probability will be calculated and so are their entropy.

Finally, we will be able to understand what is most preferred topology for synapsis, and how will this preference be influenced by the parameters of the local sysnaptosome structure (the local small rigid body).

This project also serves as general Monte Carlo simulation package for worm-like chains, with multiple additional structures, like chain thickness (volume exclusion), computational topology support (calculation of Alexander Polynomials), biased potential sampling and various sampling methods.

Written in C++, the class diagram and source code documents are available at:

http://vologodskiilab6.chem.nyu.edu/

Current version SVN checkout:

**svn checkout http://montering.googlecode.com/svn/Synapsis/branch/rev%2096_fullStaticCoreHeterogenous --[revision 109](https://code.google.com/p/montering/source/detail?r=109) montering-read-only**
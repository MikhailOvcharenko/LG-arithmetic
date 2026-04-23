# On Arithmetic Mirror Symmetry for smooth Fano fourfolds

Author: Mikhail Ovcharenko.

Supplement materials for ``On Arithmetic Mirror Symmetry for smooth Fano fourfolds''.


## Structure

The file structure is organized as follows.

1. `Magma` contains computations of the Picard--Fuchs operator and associated local monodromies from the period sequence (using the Magma code from Fanosearch project).

2. `Polymake` contains verification via Polymake of the lattice width 1 assumption on lattice Minkowski summands of Newton polytope of Laurent polynomials. It also checks that these summands are actually Minkowski-indecomposable.

3. `SageMath` contains verification of Theorem 1.2. For example, `LG_Gr25_13.ipynb` corresponds to the case of a Fano fourfold complete intersection of multidegree (1,3) in Gr(2,5). In `QuiverFlagZeroLoci.ipynb` and `SmoothToricCI.ipynb` we consider Kalasnikov's and Coates-Kasprzyk-Prince’s lists, correspondingly. `Library.sage` provides auxiliary computational tools.

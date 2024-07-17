HUBB-CODE

F90 code for Exact Diagonalization (ED) of the Hubbard Hamiltonian in Real (Slater Determinant) Space for arbitrary lattice geometries and lattice fillings with Periodic Boundary Conditions (PBC)
or Open Boundary Conditions (OBC).

- Features available for the ED in Real Space:
  - Gound State (GS) Wave Function as a CI-Expansion.
  - Spectrum of the Hamiltonian.
  - (squared) Total Spin of the GS Wave Function.
  - C-Metric of the Fermionic Sign Problem.
  - Diagonalization of the Stoquastized Hamiltonian.
  - Diagonalization of the Hamiltonian in a Subspace (truncated Hilbert Space).
  - Fixed-Node Approximation (FNA) -  only available for BULK Systems.
    - Diagonalization of the Hubbard Hamiltonian in K-Space.
    - GS Wave Function as a CI-Expansion in K-Space.
    - (squared) Total Spin of the GS Wave Function.
    - Implementation of a Multi-Configurational Trial Wave Function with desired number of K-Space Determinants and Coefficients
      obtained form the exact solution or a sub-space diagonalization.
    - Calculation of the Sign-Flip Potential of the Trial Wave Funciton and Exact Wave Funciton.
    - Relaxation of the FNA according to the P/Q Sub-Space Partition Scheme.
    - Diagonalization of the Stoquastized FNA-Hamiltonian.

- Input description:
  - SYSTEM
    - System Type                   : bulk/ladder/other
    - Boundary Conditions           : PBC[X,Y]/PBC[X]-OBC[Y]/OBC[X]-PBC[Y]/OBC[X,Y]
    - Number of Spin-Up   Electrons
    - Number of Spin-Down Electrons
    - Calculation of S^2 Matrix (activation keyword)
    - Calculation of <S^2> Spin (activation keyword)
  - LATTICE PARAMETERS
    - A1o Lattice Parameter
    - A2o Lattice Parameter
  - LATTICE VECTORS
    - A1 Lattice Vector
    - A2 Lattice Vector
  - HUBBARD PARAMETERS
    - Hopping (t)
    - On-Site Repulsion (U)
  - STOQUASTIZED
    - Diagonalization of Stoquastized Hamiltonian (activation keyword)
    - Calculation of <S^2> Spin (activation keyword)
  - HILBER SPACE
    - Truncation of Hilbert Space (activation keyword)
    - Calculation of <S^2> Spin   (activation keyword)
  - FIXED-NODE APPROXIMATION
    - FNA Calculation (activation keyword)
    - Trucation of FNA-Hamiltonian (activation keyword) - always active recommended
    - Calculation of <S^2> Spin (activation keyword)
  - TRIAL WAVE FUNCTION
    - Hopping (t)           for K-Space Hamiltonian
    - On-Site Repulsion (U) for K-Space Hamiltonian
    - Number of K-Space Determinants in the CI-Expansion of the Trial Wave Function
    - Sub-Space Diagonalization (activation keyword)
  - FNA RELAXATION
    - P/Q Partition Relaxation scheme (activation keyword)
    - Dimension of the P-Space
  - FNA STOQUASTIZED
    - Diagonalization of Stoquastized FNA-Hamiltonian (activation keyword)
    - Calculation of <S^2> Spin                       (activation keyword)

# FEM-laplacian

This repository contains a *MATLAB* function for computing the Finite Element Method (FEM) discretization of the Laplace-Beltrami operator on triangular meshes. The discretization is parametric on the polynomial's order of the hat functions.

## Requirements
The function is self contained and does not require any additional code.

## Install
The file `FEM_higher.m` with the code can be placed in any package of the project.  
The directory `Cache` must be placed in the main directory of the project.

## Signature
The function provides two different signatures:
```
    [S, A] = FEM_higher(M, p, [BoundaryCondition, [CacheDir]])
    [S, A, Al] = FEM_higher(M, p, [BoundaryCondition, [CacheDir]])
```
where
  - `M` is a triangular mesh (see [Details](#details) section);
  - `p` is the polynomial's order of the hat function;
  - `BoundaryCondition` is a string (either `"Neumann"` or `"Dirichlet"`) determining the boundary condition for partial shapes. Default value is `"Neumann"`;
  - `CacheDir` is the path to a directory containing the cached coefficients of the hat functions (they can be computed once, since they are independent on the mesh). Default value is `"Cache/FEM"`.

The output arguments are
  - `S` the stiffness matrix;
  - `A` the mass matrix;
  - `Al` the lumped mass matrix.
 

### Details
The triangular mesh `M` must be a data structure containing at least the following fields:
  - `n` the number of vertices;
  - `m` the number of triangles;
  - `VERT` a `n-by-3` matrix containing the 3D coordinates of the vertices;
  - `TRIV` a `m-by-3` matrix containing the triangulation of the mesh (each row being a triplet of indices of vertices).

The value of the stiffness and mass matrices at the real vertices are stored in consistent order with `M.VERT` in the top-left `n-by-n` sub-matrix.  

The repository already provides coefficients for hat functions up to 9th order.

# MPNA : Small linear algebra librarie

The goal of this lab is to experiment with linear solver applied on a 2D 
poisson matrix problem where the matrix is stored is CSR format. General
format was added in order to compare things out and some other formats would be 
welcomed.

## Build 

In order to build the lib, in the project directory : 
```bash
$ cmake -S . -B Build -DCMAKE_BUILD_TYPE=<type>
$ cd Build
$ make
```
and to run the test to verify everything is working flawlessly :
```bash
./main.exe
```
(! Tests needs to be fixed with new implems)

## About :
The mesh stencil for the poisson matrix is the following :
```
    -1  
  -1 4 -1  
    -1  
```
Implementend solvers :
- Jacobi
- Gauss Seidel
- Conjugate Gradient
- GMRES (relies on lapack)

Available Matrix formats :
- CSR
- General


# Ressources : 
https://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html  
https://en.wikipedia.org/wiki/Sparse_matrix  
http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf  
https://fr.wikipedia.org/wiki/Diff%C3%A9rence_finie  
https://en.wikipedia.org/wiki/Conjugate_gradient_method  
https://en.wikipedia.org/wiki/Generalized_minimal_residual_method  
https://en.wikipedia.org/wiki/Arnoldi_iteration  
https://en.wikipedia.org/wiki/Hessenberg_matrix  

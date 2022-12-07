# projetSVD
Projet de MPNA (CHPS 2022)
Implementation de la methode de décomposition en valeur singulière (SVD) d'une matrice rectanglaire représenant un dataset quelconque.

#Theorie sur la SVD
Soit une matrice réelle A de dimension $m \time n$. La SVD de A donne :  
#
hdr  main.c  makefile  obj  perf.c  proto  README.md  src  test  

- ./hdr:  
benchmark.h          gram_schmidt.h  my_seq_blas1.h  QR.h  
bidiagonalization.h  io_matrix.h     my_seq_blas2.h  SVD.h  
gen_mat.h            math_mat.h      my_seq_blas3.h  trigonalisation.h  

- ./obj:  

- ./proto:  
householder_transformation.sce  

- ./src:  
bidiagonalization.c  gram_schmidt.c  math_mat.c      my_seq_blas2.c  SVD.c  
gen_mat.c            io_matrix.c     my_seq_blas1.c  QR.c            trigonalisation.c  

- ./test:  
AtA_AAt.c  hessenberg.c  QR_decomp.c  SVD_1.c  test_mat.h  
Bidiag.c   math_mat.c    QR_method.c  SVD_3.c  test_utils.h  

# Build
- `make` : build all executables (including tests)
- `make check` : build and execute all tests
- `make clean` : remove all executables
- `make progname ; ./progname` : build and run progname.c executable

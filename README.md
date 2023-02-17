# projetSVD  
Projet de MPNA (CHPS 2022)  
Implementation de la methode de décomposition en valeur singulière (SVD) d'une matrice rectangulaire. Plusieurs algorithmes sont implémentés ici, notamment la méthode QR et l'algorithme TSQR (Demmel) pour la parallélisation, la bidiagonalisation de Golub-Kahan-Lanczos, l'utilisation des matrices de Hessenberg, et autres.
MPI est utilisé pour la parallélisation de l'algorithme TSQR.

=======
# Theorie générale
Soit une matrice réelle A de dimension $m \times n$. La SVD de A donne :  $$A=U \Sigma V^T$$  
où $U=[u_1,u_2,...,u_m] \in \mathbb{R}^{m\times n}$ et $V=[v_1,v_2,..,v_n] \in \mathbb{R}^{n\times m}$ ,sont les matrices orthogonales ( $U^T=U^{-1},V^T=V^{-1}$ ), $\Sigma$ est une matrice carré diagonale de dimension $r \times r$ tel que $\Sigma = diag(\sigma_1,\sigma_2,...,\sigma_r)$ avec $\sigma_1 \geq \sigma_2 \geq... \geq \sigma_r$ et r = min(m,n).  
En géneral les première lignes de V sont nommés vecteurs singuliers droite et gauche de A et r est le rang de la matrice A. Les $(\sigma_i)_{i \in {[1,r]}}$ sont les valeurs singulieres de A.

# Build
- `make` : build all executables (including tests)
- `make check` : build and execute all tests
- `make clean` : remove all executables
- `make progname ; ./progname` : build and run progname.c executable

# Organisation du projet
hdr  main.c  makefile  obj  perf.c  proto  README.md  src  test  

- ./hdr: Dossier des fichiers headers contenant les prototypes des fonctions et des  structures ainsi que les macros.
```   
benchmark.h          gram_schmidt.h  my_seq_blas1.h  QR.h  QRpar.h
bidiagonalization.h  io_matrix.h     my_seq_blas2.h  SVD.h  
gen_mat.h            math_mat.h      my_seq_blas3.h  trigonalisation.h  
```
- ./obj:  Dossier qui contiendra les fichiers obj et dep que la compilation produira  

- ./proto:  Dossier qui contient certains prototype des algorithme de base du projet.  
```householder_transformation.sce  ```

- ./src: Dossier contenant les fichiers sources du projet tel que le main du programme finale ainsi que les définitions des fonctions de la methode QR et de la methode SVD.
```
bidiagonalization.c  gram_schmidt.c  math_mat.c      my_seq_blas2.c  SVD.c  
gen_mat.c            io_matrix.c     my_seq_blas1.c  QR.c QRpar.h    trigonalisation.c  
```
- ./test:  Dossier contenant les tests de validations des fonctions implementer par l'equipe, il contiendra les fichiers sources ,headers et binaires des test.  
```
AtA_AAt.c  hessenberg.c  QR_decomp.c  SVD_1.c  test_mat.h  
Bidiag.c   math_mat.c    QR_method.c  SVD_3.c  test_utils.h  
```

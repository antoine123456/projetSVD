# projetSVD  
Projet de MPNA (CHPS 2022)  
Implementation de la methode de décomposition en valeur singulière (SVD) d'une matrice rectanglaire représenant un dataset quelconque.

# Theorie sur la SVD
Soit une matrice réelle A de dimension $m \times n$. La SVD de A donne :  $$A=U \Sigma V^T$$  
où $U=[u_1,u_2,...,u_m] \in \mathbb{R}^{m\times n}$ et $V=[v_1,v_2,..,v_n] \in \mathbb{R}^{n\times m}$ ,sont les matrices orthogonales ( $U^T=U^{-1},V^T=V^{-1}$ ), $\Sigma$ est une matrice carré diagonale de dimension $r \times r$ tel que $\Sigma = diag(\sigma_1,\sigma_2,...,\sigma_r)$ avec $\sigma_1 \geq \sigma_2 \geq... \geq \sigma_r$ et r = min(m,n).  
En géneral les première lignes de V sont nommés vecteurs singuliers droite et gauche de A et r est le rang de la matrice A. Les $(\sigma_i)_{i \in {[1,r]}}$ sont les valeurs singulieres de A.
# Organisation du projet
hdr  main.c  makefile  obj  perf.c  proto  README.md  src  test  

- ./hdr: Dossier des fichiers headers contenant les prototypes des fonctions et des  structures ainsi que les macros.
```   
benchmark.h          gram_schmidt.h  my_seq_blas1.h  QR.h  
bidiagonalization.h  io_matrix.h     my_seq_blas2.h  SVD.h  
gen_mat.h            math_mat.h      my_seq_blas3.h  trigonalisation.h  
```
- ./obj:  Dossier qui contiendra les fichiers obj et dep que la compilation produira  

- ./proto:  Dossier qui contient certains prototype des algorithme de base du projet.  
```householder_transformation.sce  ```

- ./src: Dossier contenant les fichiers sources du projet tel que le main du programme finale ainsi que les définitions des fonctions de la methode QR et de la methode SVD.
```
bidiagonalization.c  gram_schmidt.c  math_mat.c      my_seq_blas2.c  SVD.c  
gen_mat.c            io_matrix.c     my_seq_blas1.c  QR.c            trigonalisation.c  
```
- ./test:  Dossier contenant les tests de validations des fonctions implementer par l'equipe, il contiendra les fichiers sources ,headers et binaires des test.  
```
AtA_AAt.c  hessenberg.c  QR_decomp.c  SVD_1.c  test_mat.h  
Bidiag.c   math_mat.c    QR_method.c  SVD_3.c  test_utils.h  
```

# Alogrithmes 
## Algorithme avec la reduction de Hessenberg
D'abord on doit calculer soit les eigenvalues de la matrice $B=AA^T \in \mathbb{R}^{m\times m}$ ou $B=A^TA \mathbb{R}^{n\times n}$ selon le minimum entre m et n. En effet on a :
$$A^TA = V \Sigma U^T U \Sigma V^T = V \Sigma²V^T$$
$$AA^T = U \Sigma V^T V \Sigma U^T = U \Sigma²U^T$$
Par conséquent, U (resp.V) , $\Sigma²=\Lambda = diag(\lambda_1,...,\lambda_r)$  sont les vecteurs propres et valeurs propres de la matrice B. Autrement dit les valeurs propres de B sont les carrés des valeurs singulières de A :   $\lambda_i²_{i \in {[1,r]}}$ = $\sigma_i²_{i \in {[1,r]}}$.  
On sait que B est une matrice symétrique alors on reduisant B sous sa forme Hessenberg on obtient une matrice triangulaire supérieur et symétrique donc cela revient a une matrice tridiagonal. Puis on applique la méthode QR à la matrice tridiagonale pour calculer les valeurs propres de B ainsi l'on obtiendra les valeurs singulieres de A.
```
**Entrée :** A,matrice rectangulaire de taille (m,n)
Si m<n alors  faire :
  B=AA^T
Sinon fiare :
  B=A^TA 
On applique une Hessenberg_Reduction(B)->H
On applique une Methode_QR(H)->Sp(B)=Lambda
**Sortie :** Lambda->Sigma
```
## Algorithme avec la methode Bi-Lanczos
# Build
- `make` : build all executables (including tests)
- `make check` : build and execute all tests
- `make clean` : remove all executables
- `make progname ; ./progname` : build and run progname.c executable

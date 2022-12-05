#include"my_seq_blas1.h"

//Fonction sequentiel qui realise y = a*x+y
//Arguments:
//	un entier dim: la dimension des vecteurs x,y de type int
//	un réel alpha: coefficient a de type double
//	les vecteurs x,y: pointeur de tableau de type double
void seq_d_axpy(int dim, double alpha, double* X, double* Y){
	for(int i = 0; i < dim; i++){
		Y[i] = alpha*X[i] + Y[i];
	}
}

//Fonction sequentiel qui realise le produit scalaire
//Arguments:
//	un entier dim: la dimension des vecteurs x,y
//	des pointeurs de tableau : vecteurs x,y
//Sortie:
//	un réel result: produit scalaire de x et y
double seq_d_dotprod(int dim, double *X, double *Y){
	double result = 0.0;

	for(int i = 0 ; i<dim ; i++){
		result = X[i]*Y[i] + result;
	}
	return result;
}

//Fonction sequentiel qui realise la norme 2 de x
//Arguments:
//	un entier dim: la dimension du vecteur x
//	le pointeur de tableau de type double: vecteur x
//Sortie:
//	un réel result: norme de x
double seq_d_norm(int dim, double *X){
	return  sqrt(seq_d_dotprod(dim, X, X));
}


double seq_beta_alpha_dprod(int dim, double alpha, double *X, double *Y){
	return alpha + seq_d_dotprod(dim, X, Y);
}

//Fonction sequentiel qui realise la reduction sur vecteurs : une somme
//Arguments:
//	un entier dim: la dimension du vecteur x
//	le pointeur de tableau de type double: vecteur x
//Sortie:
//	un réel result: la somme sur les élèments x
double seq_d_reduction(int dim, double *vect){
	double sum = 0.0;
	for(int i=0; i< dim ; i++){
		sum =+vect[i];
	}
	return sum;
}

// Multiplication of a vector X by scalar alpha
double* seq_dax(int dim, double alpha,double *X){
	double *Y;
	Y  = calloc(dim, sizeof(double*));
	for(int i =0; i<dim; i++){
		Y[i] = alpha*X[i];
	}
	return Y;
}

// Set a vector of dimension dim to 0
void init_vect(int dim, double* vect){
	for(int i = 0; i<dim; i++){
		vect[i]=0.0;
	}
}

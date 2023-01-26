#include "kernel.h"
#include "jordan.h"
#include "utils.h" // TODO : replace par math_math
#include <stdlib.h>
#include <stdio.h>
/// @brief Calcul le produit matriciel entre une matrice et un vecteur, si le résultat est le vecteur nul, on est dans le noyau de la matrice
/// @param A {double *} Matrice
/// @param Vec {double *} Vecteur
/// @param n {int} taille de la matrice
/// @return {bool} true si le vecteur n'est pas dans le noyau de la matrice
bool isNotInKernel(double *A, double *Vec, int n)
{
    double *tmpRes = pdtMatVec(A, Vec, n);
    for (int i = 0; i < n; i++)
    {
        if (tmpRes[i] != 0)
        {
            free(tmpRes);
            return true;
        }
    }
    free(tmpRes);
    return false;
}

/** @brief construit une matrice carrée
 * @param Mat matrice
 * @param m dimension 1
 * @param n dimension 2
 */
void to_square_mat(double *Mat, int m, int n)
{
    if (m != n)
    {
        if (m > n)
        {
            transpose(Mat, m, n);
        }
        double *C = (double *)calloc(n * n, sizeof(double));
        for (int l = 0; l < n; l++)
        {
            for (int c = 0; c < m; c++)
            {
                C[l * m + c] = Mat[l * m + c];
            }
        }
        Mat = C;
        free(C);
    }
}
/// @brief Regarde où se trouve la valeur val dans le tableau list, en partant de l'indice watched[i] jusqu'à watched[i]+n
/// @param list liste à observer
/// @param watched valeurs déjà présentes dans la liste des valeurs propres
/// @param val valeur à regarder
/// @param n étendue à regarder
/// @return 
int listIndex(double *list, int *watched, double val, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (list[watched[i]] == val)
            return i;
    }
    return -1;
}
/// @brief Extrait des infos sur les valeurs propres
/// @param eigVal {double *} tableau des valeurs propres
/// @param n {int} nombre de valeurs propres
/// @return nEigVal {int} nombre de valeurs propres distinctes, nMult {int *} nombre de fois où la valeur propre est présente, from {int *} indice de la première occurence de la valeur propre
DiagInfo_t getMult(double *eigVal, int n)
{
    DiagInfo_t info;
    int tik = 0;
    int where = 0;
    info.mult = (int *)calloc(n, sizeof(int)); // table des mutilplicités des valeurs propres
    info.from = (int *)calloc(n, sizeof(int)); // position initiales des premieres occurences des valeurs propres
    for (int i = 0; i < n; i++)
    {
        double cur = eigVal[i]; // valeur courante

        where = listIndex(eigVal, info.from, cur, tik); // position de la valeur courante dans la liste des valeurs propres
        // printf("\n%d ", where);
        if (where > 0) // si la valeur courante est déjà présente
        {
            info.mult[where] += 1;  // on incrémente la multiplicité de la valeur courante
        }
        else
        {
            info.mult[tik] += 1; // sinon on incrémente la multiplicité de la valeur courante
            info.from[tik] = i;  // et on stocke la position de la première occurence de la valeur courante
            tik++;               // on incrémente le nombre de valeurs propres distinctes
        }
    }
    info.nEigVal = tik; // on stocke le nombre de valeurs propres distinctes
    return info;       // on retourne les infos
}
/**
 * @brief Jordanisation d'une matrice générale
 * @param Mat matrice à jordaniser
 * @param n nombre de colonnes
 * @param m nombre de lignes
 * @param V vecteur des valeurs propres
 * @return J, P matrice de jordan et matrice de passage
 */
Jordan_Bidiag_t jordan(double *Mat, int m, int n, double *V)
{
    ///initialisation///
    Jordan_Bidiag_t ret;
    double *J = (double *)malloc(n * m * sizeof(double)); // Matrice de Jordan
    double *P = (double *)malloc(n * m * sizeof(double)); // Matrice de passage
    // TODO : make square matrix
    ////////////

    // make_square_mat(Mat, m, n);

    //////extraction de valeurs ropres et infos/////
    double *eigVect = (double *)calloc(sizeof(double), n * m); // Vecteurs propres
    eigen_t eig; // mock future intruduction des propres
    eig.vectors = eigVect; // mock future intruduction des vecteurs propres
    eig.values = V; // mock future intruduction des valeurs propres
    DiagInfo_t mult = getMult(eig.values, m); // info sur les valeurs propres
    //////////// matrice pour la boucle ////////////
    double *iden = identity(n); // Matrice identité
    double *tmpKernel = (double *)malloc(n * n * sizeof(double));  // (A-val*I)^(time-1)
    double *tmpKernel2 = (double *)malloc(n * n * sizeof(double)); // (A-val*I)^(time)
    double *matMnId = (double *)malloc(n * n * sizeof(double));    // (A-val*I)
    double *tmpVec = (double *)malloc(n * sizeof(double)); // vecteur du kernel temporaire
    double *tmpRes = (double *)malloc(n * n * sizeof(double)); // pour le calcul de la matrice élevée
    double *tmpRes2 = (double *)malloc(n * sizeof(double)); // pour le calcul des vect de Jordan
    ////variables de la boucle////
    noyau_t kernel;  // ker(A-val*I)^time
    bool notInKernel = false; // si la valeur propre est dans le noyau
    int starPos = 0; // position de la valeur propre dans la matrice de Jordan
    int time = 0;    // multiplicativité de la valeur propre
    int idx = 0;     // indice de la valeur propre dans la liste "values"
    double val = 0;  // valeur propre courante
    ///// loop ///////
    for (int i = 0; i < mult.nEigVal; i++)// pour chaque valeurs propres de A
    {
        time = mult.mult[i]; // multiplicativité de la valeur propre
        idx = mult.from[i];  // indice de la valeur propre dans la liste "values"
        val = eig.values[idx]; // valeur propre courante
        starPos = time - 1 + i; // position de la valeur propre dans la matrice de Jordan
        matMnId = subMat(Mat, iden, val, n);      // (A-val*I)
        tmpKernel = powMat(matMnId, n, time - 1); // (A-val*I)^(time-1)
        tmpKernel2 = powMat(matMnId, n, time);    // (A-val*I)^(time)
        kernel = noyau(tmpKernel2, n);            // ker(A-val*I)^time
        for (int k = 0; k < kernel.size; k++) // pour chaque vecteur propre de ker(A-val*I)^time
        {
            for (int j = 0; j < n; j++) // extraction de ce vecteurs propres
            {
                tmpVec[j] = kernel.vectors[k + j * n]; // i-ème vec propres de ker(A-val*I)^time
            }
            J[starPos * n + starPos] = val; // stockage des block de jordan dans la matrice J
            if (time > 1) // évitable : si la valeur propre est multiple
            {
                notInKernel = isNotInKernel(tmpKernel, tmpVec, n); // tmpVec ∉ ker(A-Id)^(time-1) ?
                if (notInKernel) // tmpVec ∉ ker(A-Id)^(time-1)
                {
                    store_column(P, tmpVec, n, starPos, 1); // stockage des vecteur de jordan dansla matrice P
                    for (int l = starPos-time; l < starPos - 1; l++)
                    {
                        tmpRes = powMat(matMnId, n, starPos - l - 1); // (A-val*I)^(time-1)
                        tmpRes2 = pdtMatVec(tmpRes, tmpVec, n); // (A-val*I)^(time-1)*tmpVec ∉ ker(A-Id)^(time-1)
                        store_column(P, tmpRes2, n, l+1, 1); // stockage des vecteur de jordan dansla matrice P
                        J[n*(l+1)+l+1] = val; // stockage des block jordan dansla matrice J
                        J[n*(l+1)+l+2] = val; // stockage des block jordan dansla matrice J
                        }
                    break;
                }
            }
            else // si la valeur propre est simple
            {
                store_column(P, tmpVec, n, starPos, 1); // stockage des vecteur de jordan dansla matrice P
            }
        }
    }
    free(matMnId);free(tmpKernel);free(tmpKernel2);free(tmpVec );free(tmpRes );free(tmpRes2);
    //////////////////
    ret.J = J;
    ret.P = P;
    return ret;
}
void freeInfo(DiagInfo_t d)
{
    free(d.mult);
    free(d.from);
}
void freeJordan(Jordan_Bidiag_t j)
{
    free(j.J);
    free(j.P);
}
// TODO : opérations à binder avec Math
// [ ] : powMat
// [ ] : store_column
// [ ] : subMat
// [ ] : pdtMatVec
// [ ] : transpose
// [ ] : identity
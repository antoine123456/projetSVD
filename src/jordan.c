#include <jordan.h>

void to_square_mat(double *Mat, int m, int n)
{
    if (m != n)
    {
        double *C = (double *)calloc(n * n, sizeof(double));
        for (int l = 0; l < n; l++)
        {
            for (int c = 0; c < m, c++)
            {
                C[l * m + c] = Mat[l * m + c];
            }
        }
        Mat = C;
        free(C);
    }
}

int listIndex(double *list, double *watched, double val, int n)
{
    for (int i = 0; i < n, i++)
    {
        if (list[watched[i]] == val)
            return i;
    }
    return -1
}

DiagInfo_t getMult(double *eigVal, int n)
{
    DiagInfo_t info;
    int tik = 0;
    int where = 0;
    info.mult = (int *)calloc(n, sizeof(int));
    info.from = (int *)calloc(n, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        cur = eigVal[i];
        where = listIndex(eigVal, from, cur, tik);
        if (where > 0)
        {
            mult[where] += 1;
        }
        else
        {
            tik++;
            mult[tik] += 1;
            from[tik] = i;
        }
    }
    info.nEigVal = tik;
    return info;
}



void makeJP(DiagInfo_t mult, double* Mat, eigen_t eig){

}

/**
 * @brief Jordanisation d'une matrice générale
 * @param Mat matrice à jordaniser
 * @param n nombre de colonnes
 * @param m nombre de lignes
 * @return J, P matrice de jordan et matrice de passage
 */
Jordan_Bidiag_t jordan(double *Mat, int m, int n)
{
    if (m > n)
    {
        transpose(Mat);
    }
    to_square_mat(Mat);
    eigen_t eig = QR_method(Mat, n);
    DiagInfo_t mult = getMult(eig.values, n);
    makeJP(mult, Mat, eig);
}

freeInfo(DiagInfo_t i)
{
    free(d.mult);
    free(d.From);
}
freeJordan(Jordan_Bidiag_t j)
{
    free(j.J);
    free(j.P);
}
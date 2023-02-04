#pragma once

#include <QR.h>

void ThinQR(QR_t *QR, double *A, int m, int n);

double *GenData(int m, int n);

double* QRparRoot(int m, int n, int np, int rank);
void QRparOthers(int m, int n, int np, int rank);
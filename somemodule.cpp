//
// Created by frihetstegn on 5/30/23.
//

#include <iostream>

using namespace std;


void print_matrix(double **A, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) printf("%8.5f ", A[i][j]);
        printf("\n");
    }
}

void gauss(double **A, int Nf, double *c) {
    int i, j, k;
    double res;
    double tmpV;
    double **AA;
    AA = new double *[Nf];
    for (int j = 0; j < Nf; j++) AA[j] = new double[Nf + 1];
    for (int i = 0; i < Nf; i++)
        for (int j = 0; j <= Nf; j++) AA[i][j] = A[i][j];
// метод Гаусса с выбором ведущего по столбцу
    for (k = 0; k < Nf; k++) { //ищем ведущий
        int mi = k;
        double mA = 0.0;
        for (i = k; i < Nf; i++) {
            if (abs(A[i][k]) > abs(mA)) {
                mi = i;
                mA = abs(A[i][k]);
            }
        }//for i
//переставляем текущую и ведущую строки
        if (mi != k) {
            for (j = k; j <= Nf; j++) {
                tmpV = AA[k][j];
                A[k][j] = A[mi][j];
                A[mi][j] = tmpV;
            }
        }
//Гауссово исключение
        for (j = Nf; j >= k; j--) A[k][j] = A[k][j] / A[k][k];
        for (i = 0; i < Nf; i++)
            if (i != k) {
                for (j = Nf; j >= k; j--)
                    A[i][j] = A[i][j] - A[k][j] * A[i][k];
            }//if
    }//for k
//вычисляем невязку
    res = 0.0;
    for (i = 0; i < Nf; i++) {
        tmpV = 0.0;
        for (j = 0; j < Nf; j++) tmpV += AA[i][j] * A[j][Nf];
        res += abs(tmpV - AA[i][Nf]);
    }
    printf("Error of SLAE solution is %15.13f \n", res);
    for (i = 0; i < Nf; i++) c[i] = A[i][Nf];
// удаляем динамич массив из памяти
    for (int i = 0; i < Nf; i++) delete[] AA[i];
    delete[] AA;
}


void gauss2(double **A, int n, double c[]) {

    for (int k = 0; k < n; k++) {
        for (int j = k + 1; j < n; j++) {
            double dcoeff = A[j][k] / A[k][k];
            for (int i = k; i < n; i++) {
                A[j][i] -= dcoeff * A[k][i];
            }
            A[j][n] -= dcoeff * A[k][n];
        }
    }


    for (int k = 0; k < n; k++) {
        for (int j = k + 1; j < n; j++) {
            double dcoeff = A[n-j-1][n-k-1] / A[n-k-1][n-k-1];
            for (int i = k; i < n; i++) {
                A[n-j-1][n-i-1] -= dcoeff * A[n-k-1][n-i-1];
            }
            A[n-j-1][n] -= dcoeff * A[n-k-1][n];
        }
    }

    for (int i = 0; i < n; i++) c[i] = A[i][n] / A[i][i];
}


double get_step(double a, double b, int n) {
    return (b - a) / n;
}

void grid(double points[], double a, double b, int n) {
    double h = get_step(a, b, n);
    for (int i = 0; i < n + 1; i++) {
        points[i] = a + i * h;
    }
}
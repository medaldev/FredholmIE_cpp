//
// Created by frihetstegn on 5/30/23.
//
#include <iostream>
#include <tuple>
#include "somemodule.h"
using namespace std;

double K(double x1, double x2, double y1, double y2) {
    return (x1 - x2) - (y1 - y2);
}

double u0(double x1, double x2, double lam) {
    return 1 - lam * (x1 - x2);
}


//tuple<int, int> to_i1_i2 (int I, int n) {
//    int i = I % n;
//    int j = I - i * n;
//    return {i, j};
//}

int to_I(int i1, int i2, int n) {
    return i2 * n + i1;
}


double intg2(tuple<double, double> **X, int i1, int i2, int j1, int j2, double lam, int n, int nq) {
    double sum = 0.;

    auto [a11, a21] = X[i1][i2];
    auto [a12, a22] = X[i1 + 1][i2 + 1];
    double h1 = get_step(a11, a12, nq);
    double h2 = get_step(a21, a22, nq);

    auto [b11, b21] = X[j1][j2];
    auto [b12, b22] = X[j1 + 1][j2 + 1];
    double h3 = get_step(b11, b12, nq);
    double h4 = get_step(b21, b22, nq);


    for (int k1 = 0; k1 < nq; k1++) {

        double x1 = a11 + (k1 + 0.5) * h1;

        for (int k2 = 0; k2 < nq; k2++) {

            double x2 = a21 + (k2 + 0.5) * h2;

            for (int k3 = 0; k3 < nq; k3++) {

                double x3 = b11 + (k3 + 0.5) * h3;

                for (int k4 = 0; k4 < nq; k4++) {
                    double x4 = b21 + (k4 + 0.5) * h4;

                    sum += K(x1, x2, x3, x4);
                }
            }
        }
    }
    return sum * h1 * h2 * h3 * h4;

}


double intg1(tuple<double, double> **X, int i1, int i2, double lam, int n, int nq) {
    double sum = 0.;

    auto [a11, a21] = X[i1][i2];
    auto [a12, a22] = X[i1 + 1][i2 + 1];
    double h1 = get_step(a11, a12, nq);
    double h2 = get_step(a21, a22, nq);

    // cout << a11 << " " << a12 << " " << a21 << " " << a22 << "  |  \n";


    for (int k1 = 0; k1 < nq; k1++) {

        double x1 = a11 + (k1 + 0.5) * h1;

        for (int k2 = 0; k2 < nq; k2++) {

            double x2 = a21 + (k2 + 0.5) * h2;

            sum += u0(x1, x2, lam);
        }
    }
    // cout << h1 << " " << h2 << "\n";
    return sum * h1 * h2;

}



int main_lab4() {

    const int n = 10;

    const double a1 = 0., b1 = 1.;
    const double a2 = 0., b2 = 1.;


    const double lam = 0.1;
    const int nq = 10;

    int N = n * n;

    double C[N];

    auto ** X =  new tuple<double, double> * [N + 1];

    auto ** A = new double * [N];

    for(int i1=0;i1<N+1; i1++) {
        A[i1] = new double[N + 1];
        X[i1] = new tuple<double, double>[N + 1];
    }



    double h1 = get_step(a1, b1, n);
    double h2 = get_step(a2, b2, n);

    for (int i=0; i < n + 1; i++) {

        for (int j=0; j < n + 1; j++) {
            double x_1 = a1 + i * h1;
            double x_2 = b1 * j * h2;
            X[i][j] = {x_1, x_2};
        }
    }

    // cout << "ss\n";



    for(int i1=0;i1<n; i1++){
        for(int i2=0;i2<n; i2++) {

            int I = to_I(i1, i2, n);

            for(int i3=0;i3<n; i3++) {
                for (int i4 = 0; i4 < n; i4++) {
                    int J = to_I(i3, i4, n);
                    // cout << "-> " << I << " " << J << " | (" << i1 << i2 << i3 << i4 << " \n";
                    A[I][J] = (I == J) * h1 * h2 - lam * intg2(X, i1, i2, i3, i4, lam, n, nq);
                }
            }
            // cout << "---NNN---" << I << " " << N << "\n";
            A[I][N]=intg1(X, i1, i2, lam, n, nq);
        }
    }


    //print_matrix(A, N, N + 1);
    //gauss((double **)&A, n, (double *)&C);
    gauss(A, N, C);
    for (int ci = 0; ci < N; ci++) {
        cout << "C[" << ci << "] = " << C[ci] << "   ";
        delete[] A[ci];
    }


    delete[] A;

    return 0;
}



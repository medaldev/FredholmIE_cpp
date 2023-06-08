//
// Created by frihetstegn on 5/30/23.
//

#include "lab3.h"
#include <iostream>
#include <tuple>
#include "somemodule.h"
using namespace std;


double K(double s, double t) {
    return s - t;
}

double u0(double x, double lam) {
    return 1 - lam * (x - 1. / 2.);
}

double phi_0(double x, double v1, double v2) {
    if (x >= v1 && x <= v2) return (v2 - x) / (v2 - v1);
    else return 0;
}

double phi_n(double x, double v1, double v2) {
    if (x >= v1 && x <= v2) return (x - v1) / (v2 - v1);
    else return 0;
}

double phi_i(double x, double v_prev, double v_this, double v_next) {
    if ((x >= v_prev) && (x <= v_this)) {
        return (x - v_prev) / (v_this - v_prev);
    } else if ((x >= v_this) && (x <= v_next)) {
        return (v_next - x) / (v_next - v_this);
    } else {
        return 0.;
    }
}

double phi(double x, int i, double *X, int n) {
    if (i == 0) {
        return phi_0(x, X[0], X[1]);
    } else if (i == n)  {
        return phi_n(x, X[n - 1], X[n]);
    } else {
        return phi_i(x, X[i - 1], X[i], X[i + 1]);
    }
}

tuple<double, double> get_range(int i, double *X, int n) {
    if (i == 0) {
        return {X[0], X[1]};
    } else if (i == n)  {
        return {X[n - 1], X[n]};
    } else {
        return {X[i - 1], X[i + 1]};
    }
}

double intg2(double *X, int i, int j, double lam, int n, int nq) {
    double sum = 0.;

    double h1, h2;

    auto [a1, a2] = get_range(i, X, n);
    auto [b1, b2] = get_range(j, X, n);

    h1 = get_step(a1, a2, nq);
    h2 = get_step(b1, b2, nq);


    for (int i1 = 0; i1 < nq; i1++) {
        double x1 = a1 + (i1 + 0.5) * h1;
        double left = phi(x1, i, X, n) * phi(x1, j, X, n);

        for (int i2 = 0; i2 < nq; i2++) {
            double x2 = b1 + (i2 + 0.5) * h2;

            sum += -lam * K(x1, x2) * phi(x1, i, X, n) * phi(x2, j, X, n) * h2;
        }
        sum += left;
    }

    return sum * h1;
}

double intg1(double *X, int i, double lam, int n, int nq) {
    double sum = 0.;

    double h1;

    auto [a1, a2] = get_range(i, X, n);
    h1 = get_step(a1, a2, nq);

    for (int i1 = 0; i1 < nq; i1++) {
        double x1 = a1 + (i1 + 0.5) * h1;
        sum += u0(x1, lam) * phi(x1, i, X, n);
    }

    return sum * h1;
}

int main_lab3() {

    const int n = 30;
    const double a = 0., b = 1.;
    const double lam = 0.1;
    const int nq = 100;
    double *X = new double [n];

    double h = get_step(a, b, n);
    grid(X, a, b, n);



    double C[n + 1];
    double ** A = new double * [n + 1];

    for (int i = 0; i < n + 1; i++) {
        A[i] = new double [n + 2];

        for (int j = 0; j < n + 1; j++) {
            A[i][j] = intg2(X, i, j, lam, n, nq);

        }

        A[i][n + 1] = intg1(X, i, lam, n, nq);
    }
    print_matrix(A, n + 1, n + 2);
    //gauss((double **)&A, n, (double *)&C);
    gauss(A, n + 1, C);
    for (int ci = 0; ci < n + 1; ci++) {
        cout << "C[" << ci << "] = " << C[ci] << "   ";
        delete[] A[ci];
    }


    delete[] A;

    return 0;
}


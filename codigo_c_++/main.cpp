#include <iostream>
#include <cmath>

using namespace std;

void printTabela(double *dominiox, int n, double dx, double *A, double *p, double *V, double *T) {
    printf("\n\n\nRESULTADOS COMPARACAO COM TABELA 7.3 LIVRO DO Anderson\n");
    printf("____________________________________________________________________\n");
    printf("|    x/L    |    A/A*    |    p/p_0    |    V/a_0    |    T/T_0    |\n");
    printf("____________________________________________________________________\n");
    for (int i = 0; i < n; i++) {
        printf("|    %.3f  |", dominiox[0] + i * dx);
        printf("    %.3f   |", A[i]);
        printf("    %.3f    |", p[i]);
        printf("    %.3f    |", V[i]);
        printf("    %.3f    |\n", T[i]);
    }
    printf("____________________________________________________________________\n\n\n");
}

double initA(double x) {
    return 1.0 + (2.2 * ((x - 1.5) * (x - 1.5)));
}

double initP(double x) {
    return 1.0 - (0.3146 * x);
}

double initT(double x) {
    return 1.0 - (0.2314 * x);
}

double initV(double x, double T) {
    return (1.0 + (1.09 * x)) * (sqrt(T));
}

double *dpdt(double *p, double *V, double *A, int n, double dx) {
    double *result = new double[n];

    for (int i = 1; i < n - 1; i++) {
        result[i] = (-p[i] * ((V[i + 1] - V[i]) / dx))
                    - (p[i] * V[i] * ((log(A[i + 1]) - log(A[i])) / dx))
                    - (V[i] * ((p[i + 1] - p[i]) / dx));
    }

    return result;
}

double *dVdt(double *p, double *V, double *T, int n, double dx, double lambda) {
    double *result = new double[n];

    for (int i = 1; i < n - 1; i++) {
        result[i] = (-V[i] * ((V[i + 1] - V[i]) / dx))
                    - ((1 / lambda) * (((T[i + 1] - T[i]) / dx) + ((T[i] / p[i]) * ((p[i + 1] - p[i]) / dx))));
    }

    return result;
}

double *dTdt(double *V, double *T, double *A, int n, double dx, double lambda) {
    double *result = new double[n];

    for (int i = 1; i < n - 1; i++) {
        result[i] = (-V[i] * ((T[i + 1] - T[i]) / dx))
                    - (((lambda - 1) * T[i]) * (((V[i + 1] - V[i]) / dx) + (V[i] * (log(A[i + 1]) - log(A[i])) / dx)));
    }

    return result;
}

void passoPreditor(double *p, double *V, double *T, int n, double dt,
                   double *dpdt_t, double *dVdt_t, double *dTdt_t,
                   double *p_, double *V_, double *T_) {

    for (int i = 1; i < n - 1; i++) {
        p_[i] = p[i] + (dpdt_t[i] * dt);
        V_[i] = V[i] + (dVdt_t[i] * dt);
        T_[i] = T[i] + (dTdt_t[i] * dt);
    }

    p_[0] = p[0];
    p_[n - 1] = p[n - 1];
    V_[0] = V[0];
    V_[n - 1] = V[n - 1];
    T_[0] = T[0];
    T_[n - 1] = T[n - 1];
}

void
passoCorretor(double *p, double *V, double *T,
              double *A, int n, double dx, double dt, double lambda,
              double *dpdt_t, double *dVdt_t, double *dTdt_t,
              double *p_, double *V_, double *T_) {

    double *dpdt_t_dt_ = dpdt(p_, V_, A, n, dx);
    double *dVdt_t_dt_ = dVdt(p_, V_, T_, n, dx, lambda);
    double *dTdt_t_dt_ = dTdt(V_, T_, A, n, dx, lambda);

    for (int i = 1; i < n - 1; i++) {
        double dpdtAvg = 0.5 * (dpdt_t[i] + dpdt_t_dt_[i]);
        double dVdtAvg = 0.5 * (dVdt_t[i] + dVdt_t_dt_[i]);
        double dTdtAvg = 0.5 * (dTdt_t[i] + dTdt_t_dt_[i]);

        p[i] = p[i] + (dpdtAvg * dt);
        V[i] = p[i] + (dVdtAvg * dt);
        T[i] = p[i] + (dTdtAvg * dt);

    }

    delete[] dpdt_t_dt_;
    delete[] dVdt_t_dt_;
    delete[] dTdt_t_dt_;
}

void solver(double *dominiox, int n, double nt, double lambda, double courant, double dx, double dt) {

    double *A = new double[n], *p = new double[n], *V = new double[n], *T = new double[n];
    double *p_ = new double[n], *V_ = new double[n], *T_ = new double[n];

    //inicializa vetores
    for (int i = 0; i < n; i++) {

        double x = dominiox[0] + i * dx;

        A[i] = initA(x);
        p[i] = initP(x);
        T[i] = initT(x);
        V[i] = initV(x, T[i]);
    }

    double *dpdt_t, *dVdt_t, *dTdt_t;
    for (int t = 0; t < nt; t++) {

        dpdt_t = dpdt(p, V, A, n, dx);
        dVdt_t = dVdt(p, V, T, n, dx, lambda);
        dTdt_t = dTdt(V, T, A, n, dx, lambda);

        passoPreditor(p, V, T, n, dt,
                      dpdt_t, dVdt_t, dTdt_t,
                      p_, V_, T_);

        passoCorretor(p, V, T,
                      A, n, dx, dt, lambda,
                      dpdt_t, dVdt_t, dTdt_t,
                      p_, V_, T_);

        //contorno no ponto 1 (indice 0):
        p[0] = 1.0;
        T[0] = 1.0;
        V[0] = (2 * V[1] - V[2]);

        //contorno no ponto N (indice n-1):
        V[n - 1] = (2 * V[n - 2] - V[n - 3]);
        p[n - 1] = (2 * p[n - 2] - p[n - 3]);
        T[n - 1] = (2 * T[n - 2] - T[n - 3]);

        delete[] dpdt_t;
        delete[] dVdt_t;
        delete[] dTdt_t;

        printTabela(dominiox, n, dx, A, p, V, T);
    }
}


int main() {

    int n = 31;

    double dominiox[] = {0.0, 3.0};

    double courant = 0.5;

    double lambda = 1.4;

    int nt = 130;

    double dx = ((dominiox[1] - dominiox[0]) / (n - 1));

    double dt;
//    dt = courant * dx;
    dt = 2.94e-3;

    solver(dominiox, n, nt, lambda, courant, dx, dt);
    return 0;
}

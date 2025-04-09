//#define _USE_MATH_DEFINES
//
//#include <iostream>
//#include <stdio.h>
//#include <fstream>
//#include <math.h>
//#include <iomanip>
//#include <vector>
//
//using namespace std;
//
//#define eps 0.01
//
//double fi(double x, double t) {
//    return sinh(x);
//}
//
//double F(double x, double t) {
//    return 0.;
//}
//
//double psiOne(double x, double t) {
//    return exp(t);
//}
//
//double psiTwo(double x, double t) {
//    return exp(t + 1);
//}
//
//double K(double x, double t) {
//    return 1.;
//}
//
//void progonka(int N, double* A, double* B, double* C, double* f, double* y) {
//    vector<double> alpha(N + 1), beta(N + 1);
//
//    // Прямой ход
//    alpha[1] = -C[0] / B[0];
//    beta[1] = f[0] / B[0];
//
//    for (int i = 1; i < N; i++) {
//        double denominator = B[i] + A[i] * alpha[i];
//        alpha[i + 1] = -C[i] / denominator;
//        beta[i + 1] = (f[i] - A[i] * beta[i]) / denominator;
//    }
//
//    // Обратный ход
//    y[N] = (f[N] - A[N] * beta[N]) / (B[N] + A[N] * alpha[N]);
//
//    for (int i = N - 1; i >= 0; i--) {
//        y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
//    }
//}
//
//void CS(int N, double h, double tau) {
//    vector<double> A(N + 1), B(N + 1), C(N + 1), f(N + 1);
//    vector<double> x(N + 1), uPrev(N + 1), uNext(N + 1), t(N + 1);
//    double coef = 1.;
//    double ke, kw;
//
//    for (int i = 0; i < N + 1; i++) {
//        x[i] = i * h;
//        t[i] = i * tau;
//        uPrev[i] = sinh(x[i]);
//        uNext[i] = 0;
//    }
//
//    double pogr = 0;
//    float timeStart = clock() / (float)CLOCKS_PER_SEC;
//
//    for (int tt = 0; tt < N; tt++) {
//        // Граничные условия
//        A[0] = 0.;
//        B[0] = 1.;
//        C[0] = -1.;
//        f[0] = -h * exp(t[tt + 1]);
//
//        A[N] = 0.;
//        B[N] = 1.;
//        C[N] = 0.;
//        f[N] = exp(t[tt + 1] + 1);
//
//        for (int j = 1; j < N; j++) {
//            kw = (uPrev[j] + uPrev[j - 1]) / 2.;
//            ke = (uPrev[j + 1] + uPrev[j]) / 2.;
//
//            A[j] = (-coef / (h * h)) * kw;
//            B[j] = (1. / tau + coef * (ke + kw) / (h * h));
//            C[j] = (-coef / (h * h)) * ke;
//
//            f[j] = uPrev[j] / tau + (1 - coef) * (kw * uPrev[j - 1] - (ke + kw) * uPrev[j] + ke * uPrev[j + 1]) / (h * h);
//        }
//
//        progonka(N, A.data(), B.data(), C.data(), f.data(), uNext.data());
//
//        for (int i = 0; i < N + 1; i++) {
//            uPrev[i] = uNext[i];
//            double analitik = exp(t[tt + 1]) * sinh(x[i]);
//            double delta = abs(analitik - uNext[i]);
//            if (delta > pogr) {
//                pogr = delta;
//            }
//        }
//    }
//
//    float timeStop = clock() / (float)CLOCKS_PER_SEC;
//    cout << "Время = " << timeStop - timeStart << endl;
//    cout << "Ошибка = " << setprecision(10) << pogr << endl;
//}
//
//int main(int argc, char* argv[]) {
//    setlocale(LC_ALL, "RUS");
//    double h, tau;
//    double a = 0., b = 1.;
//
//    int N_values[] = { 10, 100, 1000 };
//    int num_cases = sizeof(N_values) / sizeof(N_values[0]);
//
//    for (int i = 0; i < num_cases; i++) {
//        int N = N_values[i];
//        h = (b - a) / N;
//        tau = h * h / 2;  // Условие устойчивости
//
//        cout << "N = " << N << ", h = " << h << ", t = " << tau << endl;
//        CS(N, h, tau);
//        cout << "----------------------------------------" << endl;
//    }
//
//    return 0;
//}
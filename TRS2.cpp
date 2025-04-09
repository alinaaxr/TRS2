#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include <iomanip>  // Для std::setprecision
#include <utility>  // for std::pair

using namespace std;
using namespace chrono;

// Функция для форматирования времени с высокой точностью
string format_duration(double seconds) {
    if (seconds < 1e-6) {
        return to_string(seconds * 1e9) + " мс";
    }
    //else if (seconds < 1e-3) {
    //    return to_string(seconds * 1e6) + " мс";
    //}
    //else if (seconds < 1) {
    //    return to_string(seconds * 1e3) + " мс";
    //}
    return to_string(seconds) + " с";
}

// Функция для решения трехдиагональной системы (метод прогонки)
vector<double> ThomasAlgorithm(const vector<double>& A, const vector<double>& B,
    const vector<double>& C, const vector<double>& D) {
    int n = A.size();
    vector<double> P(n), Q(n), X(n);

    P[0] = C[0] / B[0];
    Q[0] = D[0] / B[0];

    for (int i = 1; i < n; ++i) {
        double denominator = B[i] - A[i] * P[i - 1];
        P[i] = (i < n - 1) ? C[i] / denominator : 0;
        Q[i] = (D[i] - A[i] * Q[i - 1]) / denominator;
    }

    X[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        X[i] = Q[i] - P[i] * X[i + 1];
    }

    return X;
}

// Аналитическое решение
double AnalyticalSolution(double x, double t) {
    return exp(t) * sinh(x);
}

// Явная схема с расчетом ошибки и времени
pair<vector<vector<double>>, double> ExplicitScheme(double hx, double ht, int nX, int nT, double& time) {
    vector<vector<double>> u(nT, vector<double>(nX, 0));
    auto start = high_resolution_clock::now();

    // Начальное условие
    for (int k = 0; k < nX; k++) {
        u[0][k] = sinh(k * hx);
    }

    double gamma = ht / (hx * hx);

    for (int j = 1; j < nT; j++) {
        double t = j * ht;
        u[j][0] = (hx * exp(t) - u[j - 1][1]) / (hx - 1);

        for (int k = 1; k < nX - 1; k++) {
            u[j][k] = gamma * u[j - 1][k - 1] + (1 - 2 * gamma) * u[j - 1][k] + gamma * u[j - 1][k + 1];
        }

        u[j][nX - 1] = (hx * exp(t + 1) + u[j][nX - 2]) / (hx + 1);
    }

    auto stop = high_resolution_clock::now();
    time = duration_cast<microseconds>(stop - start).count() / 1e6; // в секундах

    // Расчет максимальной ошибки
    double max_error = 0;
    for (int j = 0; j < nT; j++) {
        for (int k = 0; k < nX; k++) {
            double exact = AnalyticalSolution(k * hx, j * ht);
            max_error = max(max_error, abs(u[j][k] - exact));
        }
    }

    return make_pair(u, max_error);
}

// Неявная схема с расчетом ошибки и времени
pair<vector<vector<double>>, double> ImplicitScheme(double hx, double ht, int nX, int nT, double& time) {
    vector<vector<double>> u(nT, vector<double>(nX, 0));
    auto start = high_resolution_clock::now();

    for (int k = 0; k < nX; k++) {
        u[0][k] = sinh(k * hx);
    }

    for (int j = 1; j < nT; j++) {
        double t = j * ht;
        vector<double> A(nX), B(nX), C(nX), D(nX);

        A[0] = 0; B[0] = hx - 1; C[0] = 1; D[0] = hx * exp(t);

        for (int k = 1; k < nX - 1; k++) {
            A[k] = -1.0 / (hx * hx);
            B[k] = 1.0 / ht + 2.0 / (hx * hx);
            C[k] = -1.0 / (hx * hx);
            D[k] = u[j - 1][k] / ht;
        }

        A[nX - 1] = -1; B[nX - 1] = hx + 1; C[nX - 1] = 0; D[nX - 1] = hx * exp(t + 1);

        u[j] = ThomasAlgorithm(A, B, C, D);
    }

    auto stop = high_resolution_clock::now();
    time = duration_cast<microseconds>(stop - start).count() / 1e6;

    double max_error = 0;
    for (int j = 0; j < nT; j++) {
        for (int k = 0; k < nX; k++) {
            double exact = AnalyticalSolution(k * hx, j * ht);
            max_error = max(max_error, abs(u[j][k] - exact));
        }
    }

    return make_pair(u, max_error);
}

// Схема Кранка-Николсона с расчетом ошибки и времени
pair<vector<vector<double>>, double> CrankNicolsonScheme(double hx, double ht, int nX, int nT, double& time) {
    vector<vector<double>> u(nT, vector<double>(nX, 0));
    auto start = high_resolution_clock::now();

    for (int k = 0; k < nX; k++) {
        u[0][k] = sinh(k * hx);
    }

    for (int j = 1; j < nT; j++) {
        double t = j * ht;
        vector<double> A(nX), B(nX), C(nX), D(nX);

        A[0] = 0; B[0] = hx - 1; C[0] = 1; D[0] = hx * exp(t);

        for (int k = 1; k < nX - 1; k++) {
            A[k] = -0.5 / (hx * hx);
            B[k] = 1.0 / ht + 1.0 / (hx * hx);
            C[k] = -0.5 / (hx * hx);
            D[k] = 0.5 / (hx * hx) * (u[j - 1][k - 1] + u[j - 1][k + 1]) +
                (1.0 / ht - 1.0 / (hx * hx)) * u[j - 1][k];
        }

        A[nX - 1] = -1; B[nX - 1] = hx + 1; C[nX - 1] = 0; D[nX - 1] = hx * exp(t + 1);

        u[j] = ThomasAlgorithm(A, B, C, D);
    }

    auto stop = high_resolution_clock::now();
    time = duration_cast<microseconds>(stop - start).count() / 1e6;

    double max_error = 0;
    for (int j = 0; j < nT; j++) {
        for (int k = 0; k < nX; k++) {
            double exact = AnalyticalSolution(k * hx, j * ht);
            max_error = max(max_error, abs(u[j][k] - exact));
        }
    }

    return make_pair(u, max_error);
}

// Функция для записи ошибок по пространству в фиксированный момент времени
void WriteSpaceErrors(const vector<vector<double>>& u, double hx, double ht, int fixed_t, const string& filename) {
    ofstream file(filename);
    file << "x error" << endl;
    for (int k = 0; k < u[0].size(); k++) {
        double exact = AnalyticalSolution(k * hx, fixed_t * ht);
        file << k * hx << " " << abs(u[fixed_t][k] - exact) << endl;
    }
}

// Функция для записи ошибок по времени в фиксированной точке пространства
void WriteTimeErrors(const vector<vector<double>>& u, double hx, double ht, int fixed_x, const string& filename) {
    ofstream file(filename);
    file << "t error" << endl;
    for (int j = 0; j < u.size(); j++) {
        double exact = AnalyticalSolution(fixed_x * hx, j * ht);
        file << j * ht << " " << abs(u[j][fixed_x] - exact) << endl;
    }
}

// Функция для записи норм ошибок при разных шагах
void WriteErrorNorms(const vector<double>& hx_values, const vector<double>& errors, const string& filename) {
    ofstream file(filename);
    file << "hx error" << endl;
    for (int i = 0; i < hx_values.size(); i++) {
        file << hx_values[i] << " " << errors[i] << endl;
    }
}

// Функция для записи времени выполнения
void WriteTimeResults(const vector<double>& hx_values, const vector<double>& times, const string& filename) {
    ofstream file(filename);
    file << "hx time" << endl;
    for (int i = 0; i < hx_values.size(); i++) {
        file << hx_values[i] << " " << times[i] << endl;
    }
}

int main() {
    setlocale(LC_ALL, "RUS");
    vector<double> hx_values = { 0.1, 0.01, 0.001};
    vector<double> explicit_errors, implicit_errors, cn_errors;
    vector<double> explicit_times, implicit_times, cn_times;

    cout << fixed << setprecision(6);

    // Сравнение схем для разных шагов
    for (double hx : hx_values) {
        int nX = static_cast<int>(1.0 / hx) + 1;
        double ht = hx; // Используем одинаковый шаг по времени для сравнения
        int nT = static_cast<int>(1.0 / ht) + 1;

        double time;

        cout << "\n h = " << hx << endl;

        // Явная схема
        cout << "Явная схема" << endl;
        pair<vector<vector<double>>, double> explicit_result = ExplicitScheme(hx, ht, nX, nT, time);
        explicit_errors.push_back(explicit_result.second);
        explicit_times.push_back(time);
        cout << "  Время: " << format_duration(time) <<  endl;

        // Неявная схема
        cout << "Неявная схема" << endl;
        pair<vector<vector<double>>, double> implicit_result = ImplicitScheme(hx, ht, nX, nT, time);
        implicit_errors.push_back(implicit_result.second);
        implicit_times.push_back(time);
        cout << " Время: " << format_duration(time) << endl;

        // Схема Кранка-Николсона
        cout << "Схема Кранка-Николсона" << endl;
        pair<vector<vector<double>>, double> cn_result = CrankNicolsonScheme(hx, ht, nX, nT, time);
        cn_errors.push_back(cn_result.second);
        cn_times.push_back(time);
        cout << "  Время: " << format_duration(time) << endl;

        // Запись ошибок для графиков (для последнего шага)
        if (hx == hx_values.back()) {
            WriteSpaceErrors(explicit_result.first, hx, ht, nT / 2, "explicit_space_errors.txt");
            WriteTimeErrors(explicit_result.first, hx, ht, nX / 2, "explicit_time_errors.txt");

            WriteSpaceErrors(implicit_result.first, hx, ht, nT / 2, "implicit_space_errors.txt");
            WriteTimeErrors(implicit_result.first, hx, ht, nX / 2, "implicit_time_errors.txt");

            WriteSpaceErrors(cn_result.first, hx, ht, nT / 2, "cn_space_errors.txt");
            WriteTimeErrors(cn_result.first, hx, ht, nX / 2, "cn_time_errors.txt");
        }
    }

    // Запись норм ошибок и времени выполнения
    WriteErrorNorms(hx_values, explicit_errors, "explicit_error_norms.txt");
    WriteErrorNorms(hx_values, implicit_errors, "implicit_error_norms.txt");
    WriteErrorNorms(hx_values, cn_errors, "cn_error_norms.txt");

    WriteTimeResults(hx_values, explicit_times, "explicit_times.txt");
    WriteTimeResults(hx_values, implicit_times, "implicit_times.txt");
    WriteTimeResults(hx_values, cn_times, "cn_times.txt");

    cout << "\nCalculating comparable times for target error..." << endl;
    double target_error = 0.001;
    vector<double> comparable_hx = { 0.001, 0.01, 0.1 };
    vector<double> comparable_times;

    for (double hx : comparable_hx) {
        int nX = static_cast<int>(1.0 / hx) + 1;
        double ht = hx;
        int nT = static_cast<int>(1.0 / ht) + 1;
        double time;

        cout << "Running Crank-Nicolson for hx = " << hx << "..." << endl;
        pair<vector<vector<double>>, double> cn_result = CrankNicolsonScheme(hx, ht, nX, nT, time);
        comparable_times.push_back(time);
        cout << "  Time: " << format_duration(time) << ", Max error: " << cn_result.second << endl;
    }

    WriteTimeResults(comparable_hx, comparable_times, "comparable_times.txt");

    cout << "\nAll calculations completed successfully!" << endl;
    return 0;
}

//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <string>
//#include <chrono>
//#include <iomanip>  // Для std::setprecision
//#include <utility>  // for std::pair
//
//using namespace std;
//
//using namespace chrono;
//
//// Функция для форматирования времени с высокой точностью
//string format_duration(double seconds) {
//    if (seconds < 1e-6) {
//        return to_string(seconds * 1e9) + " нс";
//    }
//    else if (seconds < 1e-3) {
//        return to_string(seconds * 1e6) + " мкс";
//    }
//    else if (seconds < 1) {
//        return to_string(seconds * 1e3) + " мс";
//    }
//    return to_string(seconds) + " с";
//}
//
//// Функция для решения трехдиагональной системы (метод прогонки)
//vector<double> ThomasAlgorithm(const vector<double>& A, const vector<double>& B,
//    const vector<double>& C, const vector<double>& D) {
//    int n = A.size();
//    vector<double> P(n), Q(n), X(n);
//
//    P[0] = C[0] / B[0];
//    Q[0] = D[0] / B[0];
//
//    for (int i = 1; i < n; ++i) {
//        double denominator = B[i] - A[i] * P[i - 1];
//        P[i] = (i < n - 1) ? C[i] / denominator : 0;
//        Q[i] = (D[i] - A[i] * Q[i - 1]) / denominator;
//    }
//
//    X[n - 1] = Q[n - 1];
//    for (int i = n - 2; i >= 0; --i) {
//        X[i] = Q[i] - P[i] * X[i + 1];
//    }
//
//    return X;
//}
//
//// Аналитическое решение
//double AnalyticalSolution(double x, double t) {
//    return exp(t) * sinh(x);
//}
//
//// Явная схема с расчетом ошибки и времени
//pair<vector<vector<double>>, double> ExplicitScheme(double hx, double ht, int nX, int nT, double& time) {
//    vector<vector<double>> u(nT, vector<double>(nX, 0));
//    auto start = high_resolution_clock::now();
//
//    // Начальное условие
//    for (int k = 0; k < nX; k++) {
//        u[0][k] = sinh(k * hx);
//    }
//
//    double gamma = ht / (hx * hx);
//
//    for (int j = 1; j < nT; j++) {
//        double t = j * ht;
//        u[j][0] = (hx * exp(t) - u[j - 1][1]) / (hx - 1);
//
//        for (int k = 1; k < nX - 1; k++) {
//            u[j][k] = gamma * u[j - 1][k - 1] + (1 - 2 * gamma) * u[j - 1][k] + gamma * u[j - 1][k + 1];
//        }
//
//        u[j][nX - 1] = (hx * exp(t + 1) + u[j][nX - 2]) / (hx + 1);
//    }
//
//    auto stop = high_resolution_clock::now();
//    time = duration_cast<microseconds>(stop - start).count() / 1e6; // в секундах
//
//    // Расчет максимальной ошибки
//    double max_error = 0;
//    for (int j = 0; j < nT; j++) {
//        for (int k = 0; k < nX; k++) {
//            double exact = AnalyticalSolution(k * hx, j * ht);
//            max_error = max(max_error, abs(u[j][k] - exact));
//        }
//    }
//
//    return make_pair(u, max_error);
//}
//
//// Неявная схема с расчетом ошибки и времени
//pair<vector<vector<double>>, double> ImplicitScheme(double hx, double ht, int nX, int nT, double& time) {
//    vector<vector<double>> u(nT, vector<double>(nX, 0));
//    auto start = high_resolution_clock::now();
//
//    for (int k = 0; k < nX; k++) {
//        u[0][k] = sinh(k * hx);
//    }
//
//    for (int j = 1; j < nT; j++) {
//        double t = j * ht;
//        vector<double> A(nX), B(nX), C(nX), D(nX);
//
//        A[0] = 0; B[0] = hx - 1; C[0] = 1; D[0] = hx * exp(t);
//
//        for (int k = 1; k < nX - 1; k++) {
//            A[k] = -1.0 / (hx * hx);
//            B[k] = 1.0 / ht + 2.0 / (hx * hx);
//            C[k] = -1.0 / (hx * hx);
//            D[k] = u[j - 1][k] / ht;
//        }
//
//        A[nX - 1] = -1; B[nX - 1] = hx + 1; C[nX - 1] = 0; D[nX - 1] = hx * exp(t + 1);
//
//        u[j] = ThomasAlgorithm(A, B, C, D);
//    }
//
//    auto stop = high_resolution_clock::now();
//    time = duration_cast<microseconds>(stop - start).count() / 1e6;
//
//    double max_error = 0;
//    for (int j = 0; j < nT; j++) {
//        for (int k = 0; k < nX; k++) {
//            double exact = AnalyticalSolution(k * hx, j * ht);
//            max_error = max(max_error, abs(u[j][k] - exact));
//        }
//    }
//
//    return make_pair(u, max_error);
//}
//
//// Схема Кранка-Николсона с расчетом ошибки и времени
//pair<vector<vector<double>>, double> CrankNicolsonScheme(double hx, double ht, int nX, int nT, double& time) {
//    vector<vector<double>> u(nT, vector<double>(nX, 0));
//    auto start = high_resolution_clock::now();
//
//    for (int k = 0; k < nX; k++) {
//        u[0][k] = sinh(k * hx);
//    }
//
//    for (int j = 1; j < nT; j++) {
//        double t = j * ht;
//        vector<double> A(nX), B(nX), C(nX), D(nX);
//
//        A[0] = 0; B[0] = hx - 1; C[0] = 1; D[0] = hx * exp(t);
//
//        for (int k = 1; k < nX - 1; k++) {
//            A[k] = -0.5 / (hx * hx);
//            B[k] = 1.0 / ht + 1.0 / (hx * hx);
//            C[k] = -0.5 / (hx * hx);
//            D[k] = 0.5 / (hx * hx) * (u[j - 1][k - 1] + u[j - 1][k + 1]) +
//                (1.0 / ht - 1.0 / (hx * hx)) * u[j - 1][k];
//        }
//
//        A[nX - 1] = -1; B[nX - 1] = hx + 1; C[nX - 1] = 0; D[nX - 1] = hx * exp(t + 1);
//
//        u[j] = ThomasAlgorithm(A, B, C, D);
//    }
//
//    auto stop = high_resolution_clock::now();
//    time = duration_cast<microseconds>(stop - start).count() / 1e6;
//
//    double max_error = 0;
//    for (int j = 0; j < nT; j++) {
//        for (int k = 0; k < nX; k++) {
//            double exact = AnalyticalSolution(k * hx, j * ht);
//            max_error = max(max_error, abs(u[j][k] - exact));
//        }
//    }
//
//    return make_pair(u, max_error);
//}
//void WriteAbsoluteDifferenceToFile(const vector<vector<double>>& numerical,
//    double hx, double ht,
//    const string& filename) {
//    ofstream file(filename);
//    if (!file.is_open()) {
//        cerr << "Error opening file: " << filename << endl;
//        return;
//    }
//
//    int nT = numerical.size();
//    int nX = numerical[0].size();
//
//    // Записываем заголовок
//    file << "# X T Absolute_Difference" << endl;
//
//    for (int j = 0; j < nT; ++j) {
//        double t = j * ht;
//        for (int k = 0; k < nX; ++k) {
//            double x = k * hx;
//            double exact = AnalyticalSolution(x, t);
//            double abs_diff = abs(exact - numerical[j][k]);  // Разность по модулю
//            file << x << " " << t << " " << abs_diff << endl;
//        }
//        file << endl; // Пустая строка для разделения временных слоев
//    }
//}
//
//int main() {
//    vector<double> hx_values = { 0.1, 0.01, 0.001 };
//
//    for (double hx : hx_values) {
//        double ht = hx;
//        int nX = static_cast<int>(1.0 / hx) + 1;
//        int nT = static_cast<int>(1.0 / ht) + 1;
//
//        double time;
//        auto cn_result = CrankNicolsonScheme(hx, ht, nX, nT, time);
//
//        string filename = "absolute_difference_hx_" + to_string(hx) + ".dat";
//        WriteAbsoluteDifferenceToFile(cn_result.first, hx, ht, filename);
//
//        cout << "Absolute difference data written to " << filename << endl;
//    }
//
//
//    return 0;
    //setlocale(LC_ALL, "RUS");
    //vector<double> hx_values = { 0.1, 0.01, 0.001 };
    //vector<double> explicit_errors, implicit_errors, cn_errors;
    //vector<double> explicit_times, implicit_times, cn_times;

    //cout << fixed << setprecision(6);

    //// Вывод заголовка таблицы
    //cout << "              " << setw(15) << "     x" << setw(15) << "t"
    //    << setw(20) << "Ошибка" << setw(20) << endl;


    //// Сравнение схем для разных шагов
    //for (double hx : hx_values) {
    //    int nX = static_cast<int>(1.0 / hx) + 1;
    //    double ht = hx; // Используем одинаковый шаг по времени для сравнения
    //    int nT = static_cast<int>(1.0 / ht) + 1;

    //    double time;

    //    // Явная схема
    //    //pair<vector<vector<double>>, double> explicit_result = ExplicitScheme(hx, ht, nX, nT, time);
    //    //explicit_errors.push_back(explicit_result.second);
    //    //explicit_times.push_back(time);

    //    //   cout << "Явная схема"
    //    //    << setw(15) << fixed << setprecision(3) << hx
    //    //    << setw(15) << fixed << setprecision(3) << ht
    //    //    << setw(20) << scientific << setprecision(3) << explicit_result.second
    //    //    << setw(20) << endl;

    //    //pair<vector<vector<double>>, double> implicit_result = ImplicitScheme(hx, ht, nX, nT, time);
    //    //implicit_errors.push_back(implicit_result.second);
    //    //implicit_times.push_back(time);
    //    //cout << "Неявная схема" << setw(15) << hx << setw(15) << ht
    //    //    << setw(20) << implicit_result.second << setw(20) << endl;

    //    // Схема Кранка-Николсона
    //    pair<vector<vector<double>>, double> cn_result = CrankNicolsonScheme(hx, ht, nX, nT, time);
    //    cn_errors.push_back(cn_result.second);
    //    cn_times.push_back(time);
    //    cout << "Явная схема" << setw(12) << hx << setw(15) << ht
    //        << setw(20) << cn_result.second << setw(20) << endl;


    //}

//    return 0;
//}
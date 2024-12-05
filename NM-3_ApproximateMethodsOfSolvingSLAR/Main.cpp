#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#define ANSI_COLOR_BLUE "\033[34m"
#define ANSI_COLOR_RESET "\033[0m"
#define ANSI_COLOR_GREEN "\033[32m"
#define ANSI_COLOR_RED "\033[31m"

using namespace std;

/**
 * @brief Функція перевірки умови збіжності для метода Якобі
 * @param matrix Матриця коефіцієнтів Х
 * @return True, якщо коефіцієнти задовольняють умові збіжностіб інакше - false
 */
bool isJacobiConvergent(const vector<vector<double>>& matrix) {

	for (int i = 0; i < matrix.size(); i++) {
		double sum = 0;
		for (int j = 0; j < matrix[i].size(); j++) {
			if (i != j) sum += abs(matrix[i][j]);
		}
		if (abs(matrix[i][i]) <= sum) return false;
	}
	return true;
}

/**
 * @brief Функція, яка розв'язує методом Якобі СЛАР 
 * @param matrix Матриця коефіцієнтів Х
 * @param B Вектор вільних членів
 * @param error Точність розрахунків
 */
void JacobiMethod(const vector<vector<double>>& matrix, const vector<double>& B, const double& error) {

	cout << "-----Welcome to Jacobi Method-----" << endl;
	if (!isJacobiConvergent(matrix)) throw "Error: <The Jacobi Method isn't convergent for this matrix>";
	int size = matrix.size();
	vector<double> Xn(size, 0);
	vector<double> Xn_new(size, 0);
	double div;
	int iteration = 1;
	cout << ANSI_COLOR_RED << setw(4) << "n" << setw(20) << "x1" << setw(20) << "x2" << setw(20) << "x3" << setw(20) << "x4" << ANSI_COLOR_RESET << endl;
	cout << setw(4) << 0 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << endl;
	do {
		for (int i = 0; i < size; i++) {
			Xn_new[i] = B[i];
			for (int j = 0; j < size; j++) {
				if (i != j) Xn_new[i] -= matrix[i][j] * Xn[j];
				else div = 1 / matrix[i][j];
			}
			Xn_new[i] *= div;
		}
		cout << setw(4) << iteration;
		for (auto i : Xn_new) {
			cout << fixed << setprecision(3) << setw(20) << i;
		}cout << endl; cout.unsetf(ios::fixed);
		double errorCurrent = 0.0;
		for (int i = 0; i < size; i++) {
			errorCurrent = max(errorCurrent, abs(Xn_new[i] - Xn[i]));
		}
		if (errorCurrent < error) return;
		iteration++;
		Xn = Xn_new;
	} while (iteration < 100);

	throw "Error: <Jacobi method didn't found Xn (iteration limits)>";
}

/**
 * @brief Функція, яка розв'язує методом Зейделя СЛАР 
 * @param matrix Матриця коефіцієнтів Х
 * @param B Вектор вільних членів
 * @param error Точність розрахунків
 */
void SeidelMethod(const vector<vector<double>>& matrix, const vector<double>& B, const double& error) {
	
	cout << "-----Welcome to Seidel Method-----" << endl;
	if (!isJacobiConvergent(matrix)) throw "Error: <The Seidel Method isn't convergent for this matrix>";
	int size = matrix.size();
	vector<double> Xn(size, 0);
	vector<double> Xn_dynamic(size, 0);
	vector<double> Xn_new(size, 0);
	double div;
	int iteration = 1;
	cout << ANSI_COLOR_RED << setw(4) << "n" << setw(20) << "x1" << setw(20) << "x2" << setw(20) << "x3" << setw(20) << "x4" << ANSI_COLOR_RESET << endl;
	cout << setw(4) << 0 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << endl;
	do {
		for (int i = 0; i < size; i++) {
			Xn_new[i] = B[i];
			for (int j = 0; j < size; j++) {
				if (i != j) Xn_new[i] -= matrix[i][j] * Xn_dynamic[j];
				else div = 1 / matrix[i][j];
			}
			Xn_new[i] *= div;
			Xn_dynamic[i] = Xn_new[i];
		}
		cout << setw(4) << iteration;
		for (auto i : Xn_new) {
			cout << fixed << setprecision(3) << setw(20) << i;
		}cout << endl; cout.unsetf(ios::fixed);
		double errorCurrent = 0.0;
		for (int i = 0; i < size; i++) {
			errorCurrent = max(errorCurrent, abs(Xn_new[i] - Xn[i]));
		}
		if (errorCurrent < error) return;
		iteration++;
		Xn = Xn_new;
	} while (iteration < 100);

	throw "Error: <Seidel method didn't found Xn (iteration limits)>";
}

int main() {

	vector<vector<double>> matrix = {
		{65, 17, -27},
		{24, -63, 27},
		{18, 23, 57}
	};
	vector<double> B = { 55, -12, 98 };

	vector<vector<double>> matrixTask = {
		{76, 21, 6, -34},
		{12, -114, 8, 9},
		{16, 24, -100, -35},
		{23, -8, 5, -75}
	};
	vector<double> BTask = { -142, 83, -121, 85 };

	try {
		JacobiMethod(matrixTask, BTask, 0.001);
	}
	catch (const char* arr) {
		cerr << arr << endl;
	}
	try {
		SeidelMethod(matrixTask, BTask, 0.001);
	}
	catch (const char* arr) {
		cerr << arr << endl;
	}


	return 0;
}
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void JacobiMethod(const vector<vector<double>>& matrix, const vector<double>& B, const double& error) {

	int size = matrix.size();
	vector<double> Xn(size, 0);
	vector<double> Xn_new(size, 0);
	double div;
	int iteration = 1;
	do {
		for (int i = 0; i < size; i++) {
			Xn_new[i] = B[i];
			for (int j = 0; j < size; j++) {
				if (i != j) Xn_new[i] -= matrix[i][j] * Xn[j];
				else div = 1 / matrix[i][j];
			}
			Xn_new[i] *= div;
		}
		cout << setw(5) << iteration;
		for (auto i : Xn_new) {
			cout << setw(20) << i;
		}cout << endl;
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

void SeidelMethod(const vector<vector<double>>& matrix, const vector<double>& B, const double& error) {
	
	int size = matrix.size();
	vector<double> Xn(size, 0);
	vector<double> Xn_dynamic(size, 0);
	vector<double> Xn_new(size, 0);
	double div;
	int iteration = 1;

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
		cout << setw(5) << iteration;
		for (auto i : Xn_new) {
			cout << setw(20) << i;
		}cout << endl;
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

int main() {

	vector<vector<double>> matrix = {
		{65, 17, -27},
		{24, -63, 27},
		{18, 23, 57}
	};
	vector<double> B = { 55, -12, 98 };

	try {
		//JacobiMethod(matrix, B, 0.001);
	}
	catch (const char* arr) {
		cerr << arr << endl;
	}
	try {
		SeidelMethod(matrix, B, 0.001);
	}
	catch (const char* arr) {
		cerr << arr << endl;
	}


	return 0;
}
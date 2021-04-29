#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;
const double PI = acos(-1);
const double tau = 2 * PI;
void FFT(vector<complex<double>> &poly) {
	// Base Case: If polynomial has only one element return
	int n = poly.size();
	if (poly.size() == 1) {
		return;
	}
	// Create two n/2 arrays for divide and conquer and store odd and even power terms seperately
	vector<complex<double>> even(n / 2), odd(n / 2);
	for (int i = 0; i < n / 2; i++) {
		even[i] = poly[i * 2];
		odd[i] = poly[(i * 2) + 1];
	}
	FFT(even);
	FFT(odd);
	// Perform the Discrete Fourier Transform operations after calculations to find DFT(even) and DFT(odd)
	double angle =  2 * PI / n;
	complex<double> w(1);
	complex<double> wn(cos(angle), sin(angle));
	for (int i = 0; i < n / 2; i++) {
		poly[i] = even[i] + (w * odd[i]);
		poly[(n / 2) + i] = even[i] - (w * odd[i]);
		w = w * wn;
	}
}
void inverseDFT(vector<complex<double>> &dftpoly) {
	// Base Case: If polynomial has only one element return
	int n = dftpoly.size();
	if (dftpoly.size() == 1) {
		return;
	}
	// Create two n/2 arrays for divide and conquer and store odd and even power terms seperately
	vector<complex<double>> even(n / 2), odd(n / 2);
	for (int i = 0; i < n / 2; i++) {
		even[i] = dftpoly[i * 2];
		odd[i] = dftpoly[(i * 2) + 1];
	}
	inverseDFT(even);
	inverseDFT(odd);
	// Perform the Inverse operations to find invDFT(even) and invDFT(odd)
	double angle = 2 * PI / n * -1;
	complex<double> w(1);
	complex<double> wn(cos(angle), sin(angle));
	for (int i = 0; i < n / 2; i++) {
		dftpoly[i] = even[i] + (w * odd[i]);
		dftpoly[(n / 2) + i] = even[i] - (w * odd[i]);
		dftpoly[i] /= 2;
		dftpoly[(n / 2) + i] /= 2;
		w = w * wn;
	}
}
vector<int> multiplication(vector<int> &A, vector<int> &B) {
	// Convert the integer coefficients into complex types
	vector<complex<double>> p1(A.begin(), A.end());
	vector<complex<double>> p2(B.begin(), B.end());
	int n = 1;
	// Resize the coefficient vector with the nearest power of two
	while (n < A.size() + B.size()) {
		n <<= 1;
	}
	p1.resize(n), p2.resize(n);
	FFT(p1);
	FFT(p2);
	for (int i = 0; i < n; i++) {
		p1[i] = p1[i] * p2[i];
	}
	inverseDFT(p1);
	vector<int> result(n);
	for (int i = 0; i < n; i++) {
		result[i] = (round(p1[i].real()));
	}
	return result;
}
int main() {
	vector<int> A, B;
	// A = 1 + 2x
	// B = 2 + 3x
	// Hence resulting polynomial C = A * B = 2 + 7x + 6x^2
	vector<int> C;
	for (int i = 0; i < 2; i++) {
		A.push_back(i + 1);
		B.push_back(i + 2);
	}
	C = multiplication(A, B);
	for (int i : C) {
		cout << i << " ";
	}
	return 0;
}

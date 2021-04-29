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
int main() {
	vector<complex<double>> a;
	for (int i = 0; i < 4; i++) {
		a.push_back(i + 1);
	}
	for (int i = 0; i < 4; i++) {
		cout << a[i] << " ";
	}
	cout << "\n";
	FFT(a);
	for (int i = 0; i < 4; i++) {
		cout << a[i] << " ";
	}
	cout << "\n";
	inverseDFT(a);
	for (int i = 0; i < 4; i++) {
		cout << a[i] << " ";
	}
	return 0;
}

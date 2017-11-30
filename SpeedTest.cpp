#include <iostream>
#include <utility>
#include <string>
#include <functional>
#include <cblas.h>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include "matrix.hpp"
#define judge

using sjtu::Matrix;

std::pair<bool, std::string> RE(const std::string &name)
{
	return std::make_pair(false, "Runtime Error! (" + name + ")");
};

std::pair<bool, std::string> WA(const std::string &name)
{
	return std::make_pair(false, "wrong answer! (" + name + ")");
};

std::pair<bool, std::string> TLE(double time)
{
	std::stringstream ss;
	ss << "Your relative time is " << time;
	return std::make_pair(false, "Time Limit Exceeded!" + ss.str());
};

std::pair<bool, std::string> PASS(double time)
{
	std::stringstream ss;
	ss << "Your relative time is " << time;
	return std::make_pair(true, ss.str());
};

void smatmul(const float *A, const float *B, float *C, int m, int p, int n, bool TA, bool TB)
{
	const CBLAS_ORDER Order = CblasRowMajor;
	const CBLAS_TRANSPOSE TransA = TA ? CblasTrans : CblasNoTrans;
	const CBLAS_TRANSPOSE TransB = TB ? CblasTrans : CblasNoTrans;
	const int M = m; //A的行数，C的行数
	const int N = n; //B的列数，C的列数
	const int K = p; //A的列数，B的行数
	const float alpha = 1;
	const float beta = 0;
	const int lda = TA ? M : K; //A的列
	const int ldb = TB ? K : N; //B的列
	const int ldc = N; //C的列
	cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void dmatmul(const double *A, const double *B, double *C, int m, int p, int n, bool TA, bool TB)
{
	const CBLAS_ORDER Order = CblasRowMajor;
	const CBLAS_TRANSPOSE TransA = TA ? CblasTrans : CblasNoTrans;
	const CBLAS_TRANSPOSE TransB = TB ? CblasTrans : CblasNoTrans;
	const int M = m; //A的行数，C的行数
	const int N = n; //B的列数，C的列数
	const int K = p; //A的列数，B的行数
	const double alpha = 1;
	const double beta = 0;
	const int lda = TA ? M : K; //A的列
	const int ldb = TB ? K : N; //B的列
	const int ldc = N; //C的列
	cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

std::pair<bool, std::string> speedMatPlus()
{
	size_t N = 1000;
	Matrix<float> a(N, N), b(N, N);
	float *A = new float[N * N], *B = new float[N * N];
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
		{
			float t = float(rand() % 10000) / 10000;
			A[i * N + j] = t;
			a(i, j) = t;
			t = float(rand() % 10000) / 10000;
			B[i * N + j] = t;
			b(i, j) = t;
		}
	int t1 = clock();
	a += b;
	int t2 = clock();
	cblas_saxpy(N * N, 1, B, 1, A, 1);
	int t3 = clock();
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			if (fabs(A[i * N + j] - a(i, j)) / A[i * N + j] > 1e-4)
			{
				std::ostringstream ss;
				ss << i << " " << j << " : " << A[i * N + j] << " " << a(i, j);
				return WA(ss.str());
			}
	double t = double(t2 - t1) / (t3 - t2);
	delete[] A;
	delete[] B;
#ifdef judge
	if (t > 3) return TLE(t);
#endif
	return PASS(t);
}

std::pair<bool, std::string> speedMul()
{
	size_t N = 1000;
	Matrix<float> a(N, N);
	float *A = new float[N * N];
	float x = float(rand() % 10000) / 100;
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
		{
			float t = float(rand() % 10000) / 10000;
			A[i * N + j] = t;
			a(i, j) = t;
		}
	int t1 = clock();
	a *= x;
	int t2 = clock();
	cblas_sscal(N * N, x, A, 1);
	int t3 = clock();
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			if (fabs(A[i * N + j] - a(i, j)) / A[i * N + j] > 1e-4)
			{
				std::ostringstream ss;
				ss << i << " " << j << " : " << A[i * N + j] << " " << a(i, j);
				return WA(ss.str());
			}
	double t = double(t2 - t1) / (t3 - t2);
	delete[] A;
#ifdef judge
	if (t > 4) return TLE(t);
#endif
	return PASS(t);
}

std::pair<bool, std::string> speedMatmul()
{
	size_t N = 300;
	Matrix<float> a(N, N), b(N, N), c(N, N);
	float *A = new float[N * N], *B = new float[N * N], *C = new float[N * N];
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
		{
			float t = float(rand() % 10000) / 10000;
			A[i * N + j] = t;
			a(i, j) = t;
			t = float(rand() % 10000) / 10000;
			B[i * N + j] = t;
			b(i, j) = t;
		}
	int t1 = clock();
	c = a * b;
	int t2 = clock();
	smatmul(A, B, C, N, N, N, false, false);
	int t3 = clock();
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			if (fabs(C[i * N + j] - c(i, j)) / C[i * N + j] > 1e-4)
			{
				std::ostringstream ss;
				ss << i << " " << j << " : " << C[i * N + j] << " " << c(i, j);
				return WA(ss.str());
			}
	double t = double(t2 - t1) / (t3 - t2);
	delete[] A;
	delete[] B;
	delete[] C;
#ifdef judge
	if (t > 40) return TLE(t);
#endif
	return PASS(t);
}


int main()
{
	std::pair<std::string, std::function<std::pair<bool, std::string>(void)>> testcases[] = {{ "speedTestMatPlus", speedMatPlus },
																							 { "speedTestMul",     speedMul },
																							 { "speedTestMatmul",  speedMatmul }};
	bool result;
	std::string information;
	for (auto &&testcase : testcases)
	{
		std::cout << testcase.first << ": ";
		std::tie(result, information) = testcase.second();
		if (result)
		{
			std::cout << "PASS. " << information << std::endl;
		} else
		{
			std::cout << "FAIL. " << information << std::endl;
		}
	}

	return 0;
}

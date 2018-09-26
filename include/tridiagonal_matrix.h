#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#define X_LEFT_BOUND 0.0
#define X_RIGHT_BOUND 1.0
#define TIME_LEFT_BOUND 0.0
#define TIME_RIGHT_BOUND 1.0

#include<vector>
class tridiagonal_matrix {
public:
	tridiagonal_matrix() {};
	~tridiagonal_matrix() {};

	tridiagonal_matrix(const unsigned int h_num, const unsigned int time_layers_num);

	long double get_h();

	std::vector<long double> TDMA();
private:
	void set_h();

	void set_tau();

	void set_a(long double alpha);

	long double f(long double x, long double t);

	long double u(long double x, long double t);

	long double mu_0(long double x);

	long double mu_1(long double t);

	long double mu_2(long double t);

	void set_above();

	void set_main();

	void set_below();

	void set_u();

	long double h = 0;
	unsigned int n = 0;
	long double tau = 0;
	long double a = 0;
	unsigned int h_num = 0;
	unsigned int time_layers_num = 0;
	std::vector<long double> above = {0};
	std::vector<long double> main = {0};
	std::vector<long double> below = {0};
	std::vector<long double> free = {0};
	std::vector<long double> u = {0};
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

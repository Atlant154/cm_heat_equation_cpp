#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#include<vector>

class tridiagonal_matrix {
public:
	/*TODO implement destructor*/
	~tridiagonal_matrix() = default;

    explicit tridiagonal_matrix(unsigned int h_num = 1, unsigned int time_layers_num = 1);

	long double get_h_();

	long double get_tau_();

	std::vector<long double> get_result_();
private:
	void set_h();

	void set_tau();

	long double f(long double x, long double t);

	long double u(long double x, long double t);

	long double mu_0(long double x);

	long double mu_1(long double t);

	long double mu_2(long double t);

	void set_above();

	void set_main();

	void set_below();

	void set_u();

	const double x_left_bound = 0.0;
	const double x_right_bound = 1.0;
	const double time_left_bound = 0.0;
	const double time_right_bound = 1.1;

	const double h = 0;
	const unsigned int n = 0;
	const double tau = 0;
	const double a = 0.0022;
	const unsigned int h_num = 0;
	const unsigned int time_layers_num = 0;
	std::vector<long double> above = {0};
	std::vector<long double> main = {0};
	std::vector<long double> below = {0};
	std::vector<long double> free = {0};
	std::vector<long double> results = {0};
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

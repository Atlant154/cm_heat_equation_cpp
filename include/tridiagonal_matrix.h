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
	double function_of_heat_sources(long double x, long double t);

	double function_of_exact_solution(long double x, long double t);

	//Time and space boundaries + diffusivity coefficient(a):
	const double x_left_bound = 0.0;
	const double x_right_bound = 1.0;
	const double time_left_bound = 0.0;
	const double time_right_bound = 1.0;
	const double a = 0.0021;

	//Unknown at compile time(define if constructor):
	const unsigned int h_num;
	const unsigned int time_layers_num;

	double h;
	double tau;
	std::vector<long double> above = {0};
	std::vector<long double> main = {0};
	std::vector<long double> below = {0};
	std::vector<long double> free = {0};
	std::vector<long double> results = {0};
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

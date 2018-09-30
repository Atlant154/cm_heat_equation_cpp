#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#include<vector>

class tridiagonal_matrix {
public:
	/*TODO implement destructor*/
	~tridiagonal_matrix() = default;

	explicit tridiagonal_matrix(unsigned int h_num = 1, unsigned int time_layers_num = 1);

	long double get_h_() const;

	long double get_tau_() const;

	std::vector<double> get_result_();
private:
	double function_of_heat_sources(double x, double t);

	double function_of_exact_solution(double x, double t);

    void get_local_result_();

	//Time and space boundaries + diffusivity coefficient(a):
	const double x_left_bound = 0.0;
	const double x_right_bound = 1.0;
	const double time_left_bound = 0.0;
	const double time_right_bound = 1.0;
	const double a = 0.0020;

	//Unknown at compile time(define in constructor):
	const unsigned int h_num;
	const unsigned int time_layers_num;


	double h;
	double tau;
	double above_coefficient;
	double main_coefficient;
	double below_coefficient;
	std::vector<double> free_part;
	std::vector<double > local_result;
	std::vector<double> results;
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

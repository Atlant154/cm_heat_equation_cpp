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
	const double x_left_bound_ = 0.0;
	const double x_right_bound_ = 1.0;
	const double time_left_bound_ = 0.0;
	const double time_right_bound_ = 1.0;
	const double a_ = 0.0020;

	//Unknown at compile time(define in constructor):
	const unsigned int h_num_;
	const unsigned int time_layers_num_;


	double h_;
	double tau_;
	double above_coefficient_;
	double main_coefficient_;
	double below_coefficient_;
	std::vector<double> free_part_;
	std::vector<double > local_result_;
	std::vector<double> results_;
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

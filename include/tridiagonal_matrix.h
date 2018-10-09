#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#include <cmath>
#include <vector>

class tridiagonal_matrix {
public:

	explicit tridiagonal_matrix(unsigned int h_num = 1, unsigned int time_layers_num = 1);

	void write_result() const;

	double get_max_error() const;

    void get_result_();

    double get_error();

private:

	double function_of_heat_sources(double x, double t);

	double function_of_exact_solution(double x, double t);

    std::vector<double> get_time_layer_result_();

	//Time and space boundaries + diffusivity coefficient(a):
	const double x_left_bound_ = 0.0;
	const double x_right_bound_ = 1.0;
	const double time_left_bound_ = 0.0;
	const double time_right_bound_ = 1.0;
	const double a_ = 0.020;

	//Unknown at compile time(define in constructor):
	const unsigned int h_num_;
	const unsigned int time_layers_num_;

	//Get the values during the initialization of
	//h_hum and time_layers_num.
	const double h_;
    const double tau_;
    const double matrix_above_;
	const double matrix_main_;

	double error_ = 0.0;
	std::vector<double> free_part_;
	std::vector<double> local_result_;
	std::vector<double> ta_result;
	std::vector<std::vector<double>> results_;
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H

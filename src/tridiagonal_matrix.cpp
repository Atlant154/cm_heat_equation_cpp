#include "../include/tridiagonal_matrix.h"

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num) : h_num(h_num),
																						   time_layers_num(time_layers_num)
{
	this -> h = (x_right_bound - x_left_bound) / h_num;
	this -> tau = (time_right_bound - time_left_bound) / time_layers_num;
	above.reserve(h_num - 1);
	main.reserve(h_num);
	below.reserve(h_num - 1);
}

long double tridiagonal_matrix::get_h_() const {
	return h;
}

long double tridiagonal_matrix::get_tau_() const {
    return tau;
}

double tridiagonal_matrix::function_of_heat_sources(long double x, long double t)
{
	return 1;
}

double tridiagonal_matrix::function_of_exact_solution(long double x, long double t)
{
	return 1;
}

std::vector<long double> tridiagonal_matrix::get_result_()
{
	std::vector<long double> alpha;
	alpha.reserve(h_num-1);
	std::vector<long double> beta;
	beta.reserve(h_num);
	std::vector<long double> local_result;
	local_result.reserve(h_num);
	alpha[0] = below[0] / main[0];
	beta[0] = free[0] / main[0];
	//Get method coefficient:
	for(auto i = 1; i < h_num - 1; ++i) {
		alpha[i] = below[i] / (main[i] - above[i - 1] * alpha[i - 1]);
		beta[i] = (free[i] - above[i - 1] * beta[i - 1]) / (main[i] - above[i - 1] * alpha[i - 1]);
	}
	beta[h_num - 1] = (free[h_num - 1] - above[h_num - 2] * beta[h_num - 2]) / (main[h_num - 1] - above[h_num - 2] * alpha[h_num - 2]);
	//Find result on time layer:
	local_result[h_num - 1] = beta[h_num - 1];
	for(auto j = h_num - 1; j > 0; --j) {
		local_result[j - 1] = beta[j - 1] - alpha[j - 1] * local_result[j];
	}
	return local_result;
}

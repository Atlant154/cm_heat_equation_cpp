#include "../include/tridiagonal_matrix.h"
#include <cmath>
#include <iostream>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num) : h_num(h_num),
																						   time_layers_num(time_layers_num)
{
	h = (x_right_bound - x_left_bound) / h_num;
	tau = (time_right_bound - time_left_bound) / time_layers_num;
	below_coefficient = above_coefficient = ((-a) * tau) / (h * h);
	main_coefficient = 1 + (2 * a * tau) / (h * h);
	free_part.reserve(h_num + 1);
	results.reserve((h_num) * (time_layers_num));
	local_result.reserve(h_num);
}

long double tridiagonal_matrix::get_h_() const {
	return h;
}

long double tridiagonal_matrix::get_tau_() const {
    return tau;
}

double tridiagonal_matrix::function_of_heat_sources(double x, double t)
{
	return 2 * t - exp(x) + x - a * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double tridiagonal_matrix::function_of_exact_solution(double x, double t)
{
	return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

std::vector<double> tridiagonal_matrix::get_result_()
{
    //Obtain solution at the initial time:
    for(auto iter = 0; iter < h_num + 1; ++iter) {
		local_result.emplace_back(function_of_exact_solution(iter * h, 0.0));
		std::cout << "X = " << iter * h << ". Exact solution = " << function_of_exact_solution(iter * h, 0.0) << ".\n";
	}
	results.insert(results.end(), local_result.begin(), local_result.end());
    for(auto iter = local_result.begin(); iter != local_result.end(); ++iter)
    	std::cout << "Local result: " << *iter << std::endl;
    for(auto time_iter = 1; time_iter < time_layers_num + 1; ++time_iter)
	{
    	for(auto space_iter = 0; space_iter < h_num + 1; ++space_iter)
    		free_part.emplace_back(local_result[space_iter] + tau * function_of_heat_sources(space_iter * h, time_iter * tau));
		//Get boundary condition:
    	free_part[h_num] += function_of_exact_solution(x_right_bound, time_iter * tau);
    	free_part[0] += function_of_exact_solution(x_left_bound, time_iter * tau);
    	get_local_result_();
		results.insert(std::end(results), std::begin(local_result), std::end(local_result));
	}
	return results;
}

void tridiagonal_matrix::get_local_result_()
{
	std::vector<double> alpha;
	alpha.reserve(h_num-1);
	std::vector<double> beta;
	beta.reserve(h_num);
	alpha.emplace_back(below_coefficient / main_coefficient);
	beta.emplace_back(free_part[0] / main_coefficient);
	//Get method coefficient:
	for(auto i = 1; i < h_num - 1; ++i) {
		alpha.emplace_back(below_coefficient / (main_coefficient - above_coefficient * alpha[i - 1]));
		beta.emplace_back((free_part[i] - above_coefficient * beta[i - 1]) / (main_coefficient - above_coefficient * alpha[i - 1]));
	}
	beta[h_num - 1] = (free_part[h_num - 1] - above_coefficient * beta[h_num - 2]) / (main_coefficient - above_coefficient * alpha[h_num - 2]);
	//Find result on time layer:
	local_result[h_num - 1] = beta[h_num - 1];
	for(auto j = h_num - 1; j > 0; --j) {
		local_result.emplace_back(beta[j - 1] - alpha[j - 1] * local_result[j]);
	}
}

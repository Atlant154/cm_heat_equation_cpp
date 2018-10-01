#include "../include/tridiagonal_matrix.h"
#include <cmath>
#include <iostream>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
: h_num_(h_num),
time_layers_num_(time_layers_num),
h_((x_right_bound_ - x_left_bound_) / h_num),
tau_((time_right_bound_ - time_left_bound_) / time_layers_num),
above_coefficient_((-a_) * tau_ / (h_ * h_)),
below_coefficient_(above_coefficient_),
free_part_(h_num_),
results_(h_num_ * time_layers_num_),
local_result_(h_num_)
{
}

long double tridiagonal_matrix::get_h_() const {
	return h_;
}

long double tridiagonal_matrix::get_tau_() const {
    return tau_;
}

double tridiagonal_matrix::function_of_heat_sources(double x, double t)
{
	return 2 * t - exp(x) + x - a_ * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double tridiagonal_matrix::function_of_exact_solution(double x, double t)
{
	return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

std::vector<double> tridiagonal_matrix::get_result_()
{
    //Obtain solution at the initial time:
    for(auto iter = 0; iter < h_num_; ++iter) {
		local_result_.emplace_back(function_of_exact_solution(iter * h_, 0.0));
		std::cout << "X = " << iter * h_ << ". Exact solution = " << function_of_exact_solution(iter * h_, 0.0) << ".\n";
	}
	results_.insert(results_.end(), local_result_.begin(), local_result_.end());
    for(auto iter = local_result_.begin(); iter != local_result_.end(); ++iter)
    	std::cout << "Local result: " << *iter << std::endl;
    for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
	{
    	for(auto space_iter = 0; space_iter < h_num_; ++space_iter)
    		free_part_.emplace_back(local_result_[space_iter] + tau_ * function_of_heat_sources(space_iter * h_, time_iter * tau_));
		//Get boundary condition:
    	free_part_[h_num_] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_right_bound_, time_iter * tau_);
    	free_part_[0] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_left_bound_, time_iter * tau_);
    	get_local_result_();
		results_.insert(std::end(results_), std::begin(local_result_), std::end(local_result_));
	}
	return results_;
}

void tridiagonal_matrix::get_local_result_()
{
	std::vector<double> alpha(h_num_);
	std::vector<double> beta(h_num_);
	alpha[0] = below_coefficient_ / main_coefficient_;
	beta[0] = free_part_[0] / main_coefficient_;
	//Get method coefficient:
	for(auto iter = 1; iter < h_num_ - 1; ++iter) {
		alpha[iter] = below_coefficient_ / (main_coefficient_ - above_coefficient_ * alpha[iter - 1]);
		beta[iter] = (free_part_[iter] - above_coefficient_ * beta[iter - 1]) / (main_coefficient_ - above_coefficient_ * alpha[iter - 1]);
	}
	beta[h_num_ - 1] = (free_part_[h_num_ - 1] - above_coefficient_ * beta[h_num_ - 2]) / (main_coefficient_ - above_coefficient_ * alpha[h_num_ - 2]);
	//Find result on time layer:
	local_result_[h_num_ - 1] = beta[h_num_ - 1];
	for(auto iter = h_num_ - 1; iter > 0; --iter) {
		local_result_[iter] = (beta[iter - 1] - alpha[iter - 1] * local_result_[iter]);
	}
}

#include "../include/tridiagonal_matrix.h"
#include <cmath>
#include <iostream>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
: h_num_(h_num + 1),
time_layers_num_(time_layers_num + 1),
h_((x_right_bound_ - x_left_bound_) / h_num),
tau_((time_right_bound_ - time_left_bound_) / time_layers_num),
above_coefficient_((-a_) * tau_ / (h_ * h_)),
below_coefficient_(above_coefficient_)
{}

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

	results_.reserve(h_num_ * time_layers_num_);
	free_part_.reserve(h_num_);
    //Obtain solution at the initial time:
    for(auto iter = 0; iter < h_num_; ++iter)
    	local_result_.emplace_back(function_of_exact_solution(iter * h_, 0.0));
	results_.insert(results_.end(), local_result_.begin(), local_result_.end());
	//Obtain solution at the remaining time:
    for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
	{
    	free_part_.clear();
    	//Get free part of equation:
    	for(auto space_iter = 0; space_iter < h_num_; ++space_iter)
    		free_part_.emplace_back(local_result_[space_iter] + tau_ * function_of_heat_sources(space_iter * h_, time_iter * tau_));
		//Get boundary condition:
    	free_part_[h_num_ - 1] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_right_bound_, time_iter * tau_);
    	free_part_[0] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_left_bound_, time_iter * tau_);
    	get_local_result_();
    	std::cout << "Time: " << time_iter * tau_ << ".\n";
    	for(auto iter = local_result_.begin(); iter != local_result_.end(); ++iter)
    		std::cout << "Value: " << * iter << ".\n";
		results_.insert(std::end(results_), std::begin(local_result_), std::end(local_result_));
	}
	return results_;
}

void tridiagonal_matrix::get_local_result_()
{
    //Get rank of free part of equation:
    auto bound = free_part_.size();
    //Declare method coefficient vectors:
	std::vector<double> alpha;
	alpha.reserve(bound);
	std::vector<double> beta(h_num_);
	beta.reserve(bound);
	alpha.emplace_back(below_coefficient_ / main_coefficient_);
	beta.emplace_back(free_part_[0] / main_coefficient_);
	//Get method coefficient:
	for(auto iter = 1; iter < h_num_ - 1; ++iter)
	{
		alpha.emplace_back(above_coefficient_ / (- main_coefficient_ - below_coefficient_ * alpha[iter - 1]));
		beta.emplace_back((below_coefficient_ * beta[iter - 1] - free_part_[iter]) / (- main_coefficient_ - below_coefficient_ * alpha[iter - 1]));
	}
	local_result_.clear();
	local_result_.reserve(bound);
	//Find result on time layer:
	local_result_.emplace_back((below_coefficient_ * beta[bound - 2] - free_part_[bound - 1]) / (- main_coefficient_ - below_coefficient_ * alpha[bound - 2]));
	for(auto iter = 1; iter < bound; ++iter) {
		local_result_.emplace_back(alpha[bound - iter - 1] * local_result_[iter - 1] + beta[bound - iter - 1]);
	}
}

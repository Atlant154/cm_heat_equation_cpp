#include "../include/tridiagonal_matrix.h"
#include <cmath>
#include <iostream>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
: h_num_(h_num + 1),
time_layers_num_(time_layers_num + 1),
h_((x_right_bound_ - x_left_bound_) / h_num),
tau_((time_right_bound_ - time_left_bound_) / time_layers_num),
above_coefficient_((-a_) * tau_ / (h_ * h_)),
below_coefficient_(above_coefficient_),
main_coefficient_(1 + (2 * a_ * tau_) / (h_ * h_)),
local_result_(h_num_),
free_part_(h_num_)
{}

long double tridiagonal_matrix::get_h() const {
	return h_;
}

long double tridiagonal_matrix::get_tau() const {
    return tau_;
}

double tridiagonal_matrix::function_of_heat_sources(double x, double t) {
    return 2 * t - exp(x) + x - a_ * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double tridiagonal_matrix::function_of_exact_solution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void tridiagonal_matrix::get_result_()
{
	for(auto iter = 0; iter < h_num_; ++iter)
		local_result_[iter] = function_of_exact_solution(iter * h_, 0.0);
	results_.emplace_back(local_result_);

	for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
	{
		for(auto space_iter = 0; space_iter < h_num_; ++space_iter)
			free_part_[space_iter] = local_result_[space_iter] + tau_ * function_of_heat_sources(space_iter * h_, time_iter * tau_);
		free_part_[0] += ((tau_ * a_) / (h_ * h_)) * function_of_exact_solution(x_left_bound_, time_iter * tau_);
		free_part_[h_num_] += ((tau_ * a_) / (h_ * h_)) * function_of_exact_solution(x_right_bound_, time_iter * tau_);
		get_time_layer_result_();
		results_.emplace_back(local_result_);
	}
}

void tridiagonal_matrix::get_time_layer_result_()
{
	auto bound = free_part_.size();
	std::vector<double> alpha(bound), beta(bound);

	alpha[0] = above_coefficient_ / main_coefficient_;
	beta[0] = free_part_[0] / main_coefficient_;

	double common_factor;

	for(auto iter = 1; iter < bound - 1; ++iter)
    {
	    common_factor = 1.0 / (main_coefficient_ - below_coefficient_ * alpha[iter - 1]);
	    alpha[iter] = above_coefficient_ * common_factor;
	    beta[iter] = (free_part_[iter] - below_coefficient_ * beta[iter - 1]) * common_factor;
    }
    beta[bound - 1] = (free_part_[bound - 1] - below_coefficient_ * beta[bound - 2])
            / (main_coefficient_ - below_coefficient_ * alpha[bound - 2]);
	local_result_[bound - 1]  = beta[bound - 1];
	for(auto iter = bound - 1; iter > 0; --iter)
	    local_result_[iter - 1] = beta[iter - 1] - alpha[iter - 1] * local_result_[iter];
}

void tridiagonal_matrix::write_result() const
{
	std::fstream result_file;
	result_file.open("../result/result.txt", std::ios::out | std::ios::trunc);
	result_file << "[[" << x_right_bound_ << "," << x_left_bound_ << ", " << h_num_ << "],"
	<< "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],";
	for(auto time_iter = 0; time_iter < time_layers_num_; ++time_iter)
	{
		result_file << "[";
		for(auto space_iter = 0; space_iter < h_num_; ++space_iter)
		{
			result_file << results_[time_iter][space_iter];
			if(space_iter + 1 != h_num_)
				result_file << ",";
		}
		result_file << "]";
		if(time_iter + 1 != time_layers_num_)
			result_file << ",";
	}
	result_file << "]\n";
	result_file.close();
}
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
	for(auto iter = 0; iter < h_num_; ++iter)
		local_result_.emplace_back(function_of_exact_solution(iter * h_num_, 0.0));

	for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
	{
		for(auto space_iter = 0; space_iter < h_num_; ++space_iter)
			free_part_.emplace_back(local_result_[space_iter] + tau_ + function_of_heat_sources(space_iter * h_, time_iter * tau_));
		free_part_[0] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_left_bound_, time_iter * tau_);
		free_part_[h_num_ - 1] += (tau_ * a_ / (h_ * h_)) * function_of_exact_solution(x_right_bound_, time_iter * tau_);
		get_local_result_();
		results_.insert(results_.end(), local_result_.begin(), local_result_.end());
		local_result_.clear();
	}
	return results_;
}

void tridiagonal_matrix::get_local_result_() {
	auto bound = local_result_.size();
	std::vector<double> alpha, beta;
	alpha.reserve(bound);
	beta.reserve(bound);

	std::cout << "Below = " << below_coefficient_ << ". Main = " << main_coefficient_ << ". Above = " << above_coefficient_ << ". \n";

	alpha.push_back(below_coefficient_ / -main_coefficient_);
	beta.push_back(free_part_[0] / main_coefficient_);

	double common_multiplyer;

	for (auto iter = 1; iter < bound - 1; ++iter)
	{
		common_multiplyer = 1.0 / ((-1.0) * main_coefficient_ - below_coefficient_ * alpha[iter - 1]);
		alpha.push_back(above_coefficient_ * common_multiplyer);
		beta.push_back((below_coefficient_ * beta[iter - 1] - free_part_[iter]) * common_multiplyer);
		std::cout << "Free part = " << free_part_[iter] << ". \n";
		std::cout << "Alpha = " << alpha[iter] << ". Beta = " << beta[iter] << ". \n";
	}

	local_result_.clear();
	local_result_.reserve(bound);

	local_result_.push_back((below_coefficient_ * beta[bound - 2] - free_part_[bound - 1]) / ((-1.0) * main_coefficient_ - below_coefficient_ * alpha[bound - 2]));

	for(auto iter = 1; iter < bound; ++iter)
		local_result_.push_back(alpha[bound - iter - 1] * local_result_[iter - 1] + beta[bound - iter - 1]);
}

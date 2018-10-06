#include "../include/tridiagonal_matrix.h"
#include <cmath>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
: h_num_(h_num + 1),
time_layers_num_(time_layers_num + 1),
h_((x_right_bound_ - x_left_bound_) / h_num),
courant_number_((-a_) / std::pow(h_, 2)),
tau_((time_right_bound_ - time_left_bound_) / time_layers_num),
above_coefficient_(tau_ * courant_number_),
below_coefficient_(above_coefficient_),
main_coefficient_(1 + 2 * courant_number_ * tau_)
{}

double tridiagonal_matrix::function_of_heat_sources(double x, double t) {
    return 2 * t - exp(x) + x - a_ * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double tridiagonal_matrix::function_of_exact_solution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void tridiagonal_matrix::get_result_()
{
    local_result_.reserve(h_num_);
    free_part_.reserve(h_num_ - 2);
	for(auto iter = 0; iter < h_num_; ++iter)
		local_result_.emplace_back(function_of_exact_solution(iter * h_, 0.0));
	results_.emplace_back(local_result_);

	for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
	{
		for(auto space_iter = 0; space_iter < h_num_ - 2; ++space_iter)
			free_part_.emplace_back(local_result_[space_iter + 1] + tau_ * function_of_heat_sources((space_iter + 1) * h_, time_iter * tau_));
		free_part_[0] += tau_ * courant_number_ * function_of_exact_solution(x_left_bound_, time_iter * tau_);
		free_part_.back() += tau_ * courant_number_ * function_of_exact_solution(x_right_bound_, time_iter * tau_);
		local_result_.clear();
		local_result_.emplace_back(function_of_exact_solution(x_right_bound_, time_iter * tau_));
		ta_result = get_time_layer_result_();
		free_part_.clear();
		local_result_.insert(local_result_.end(), ta_result.begin(), ta_result.end());
        local_result_.emplace_back(function_of_exact_solution(x_left_bound_, time_iter * tau_));
		results_.emplace_back(local_result_);
	}
}

std::vector<double> tridiagonal_matrix::get_time_layer_result_()
{
    auto bound = free_part_.size() - 1;

    std::vector<double> alpha, beta;

    alpha.reserve(bound);
    beta.reserve(bound + 1);

    alpha.emplace_back(below_coefficient_ / -main_coefficient_);
    beta.emplace_back(free_part_[0] / main_coefficient_);

    double common_factor;

    for(auto iter = 1; iter < bound - 1; ++iter)
    {
        common_factor = 1.0 / (-main_coefficient_ - below_coefficient_ * alpha.back());
        alpha.emplace_back(above_coefficient_ * common_factor);
        beta.emplace_back((below_coefficient_ * beta.back() - free_part_[iter]) * common_factor);
    }

    std::vector<double> local_solution;

    local_solution.reserve(bound + 1);

    local_solution.emplace_back((below_coefficient_ * beta.back() - free_part_.back())
    / (-main_coefficient_ - below_coefficient_ * alpha.back()));

    for(auto iter = 1; iter < bound; ++iter)
        local_solution.emplace_back(alpha[bound - iter - 1] * local_solution.back() + beta[bound - iter - 1]);

    return local_solution;
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

double tridiagonal_matrix::get_error()
{
    double abs_exact_solution, abs_our_solution, results_difference;
    for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
        for(auto space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
        {
            abs_exact_solution = std::abs(function_of_exact_solution(space_iter * h_, time_iter * tau_));
            abs_our_solution = std::abs(results_[time_iter][space_iter]);
            results_difference = std::abs(abs_exact_solution - abs_our_solution);
            if(results_difference > error_)
                error_ = results_difference;
        }
    return error_;
}

double tridiagonal_matrix::get_max_error() const
{
    return tau_ + h_ * h_;
}
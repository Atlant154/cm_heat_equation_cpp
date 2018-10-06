#include "../include/tridiagonal_matrix.h"
#include <cmath>
#include <algorithm>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
:
h_num_(h_num + 1),
time_layers_num_(time_layers_num + 1),
h_((x_right_bound_ - x_left_bound_) / h_num_),
tau_((time_right_bound_ - time_left_bound_) / time_layers_num_),
courant_number_((-a_) / std::pow(h_, 2)),
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
    local_result_.clear();
    local_result_.reserve(h_num_);
    results_.reserve(h_num_ * time_layers_num_);

    for(auto iter = 0; iter < h_num_; ++iter)
        local_result_.emplace_back(function_of_exact_solution(iter * h_, 0.0));

    results_.emplace_back(local_result_);
    free_part_.clear();

    for(auto time_iter = 1; time_iter < time_layers_num_; ++time_iter)
    {
        //Get free part of equation:
        for(auto space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            free_part_.emplace_back(local_result_[space_iter] + tau_ * function_of_heat_sources(space_iter * h_, time_iter * tau_));

        //Add bound conditions:
        free_part_[0] += courant_number_ * tau_ * function_of_exact_solution(x_left_bound_, time_iter * tau_);
        free_part_[h_num_ - 2] += courant_number_ * tau_ * function_of_exact_solution(x_right_bound_, time_iter * tau_);

        //Clear local(time layer) result vector:
        local_result_.clear();

        //Get part of result from modified Thomas Algorithm:
        ta_result = get_time_layer_result_();

        //Clear free part vector:
        free_part_.clear();

        local_result_.emplace_back(function_of_exact_solution(x_left_bound_, time_iter * tau_));
        local_result_.insert(local_result_.end(), ta_result.begin(), ta_result.end());
        local_result_.emplace_back(function_of_exact_solution(x_right_bound_, time_iter * tau_));

        results_.emplace_back(local_result_);
    }
    local_result_.clear();
}

std::vector<double> tridiagonal_matrix::get_time_layer_result_()
{
    auto n = free_part_.size();

    std::vector<double> alpha(n - 1), beta(n);
    std::vector<double> result(n);

    alpha[0] = above_coefficient_ / main_coefficient_;
    beta[0] = free_part_[0] / main_coefficient_;

    double co_factor;

    for(auto iter = 1; iter < n - 1; ++iter)
    {
        co_factor = 1.0 / (main_coefficient_ - below_coefficient_ * alpha[iter - 1]);
        alpha[iter] = above_coefficient_ * co_factor;
        beta[iter] = (free_part_[iter] - above_coefficient_ * beta[iter - 1]) * co_factor;
    }

    result[n - 1] = (free_part_[n - 1] - above_coefficient_ * beta[n - 2]) / (main_coefficient_ - below_coefficient_ * alpha[n - 2]);

    for(auto iter = n - 1; iter > 0; --iter)
        result[iter - 1] = beta[iter - 1] - alpha[iter - 1] * result[iter];

    return result;
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
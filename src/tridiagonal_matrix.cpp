#include "../include/tridiagonal_matrix.h"

#include <fstream>

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num)
        :
        h_num_(h_num + 1),
        time_layers_num_(time_layers_num + 1),
        h_((x_right_bound_ - x_left_bound_) / h_num_),
        tau_((time_right_bound_ - time_left_bound_) / time_layers_num_),
        matrix_above_((-a_ * tau_) / pow(h_, 2)),
        matrix_main_(1.0 + 2.0 * a_ * tau_ / pow(h_, 2)) {}

double tridiagonal_matrix::function_of_heat_sources(double x, double t) {
    return 2 * t - exp(x) + x - a_ * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double tridiagonal_matrix::function_of_exact_solution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void tridiagonal_matrix::get_result() {
    local_result_.reserve(h_num_);
    results_.reserve(h_num_ * time_layers_num_);

    for (unsigned iter = 0; iter < h_num_; ++iter)
        local_result_.emplace_back(function_of_exact_solution(iter * h_, 0.0));

    results_.emplace_back(local_result_);
    free_part_.clear();

    const double courant_number = (a_ * tau_) / pow(h_, 2);

    for (unsigned time_iter = 1; time_iter <= time_layers_num_; ++time_iter) {
        for (unsigned space_iter = 1; space_iter < h_num_; ++space_iter)
            free_part_.emplace_back(
                    local_result_[space_iter] + tau_ * function_of_heat_sources(space_iter * h_, time_iter * tau_));

        free_part_.front() += courant_number * function_of_exact_solution(x_left_bound_, time_iter * tau_);
        free_part_.back() += courant_number * function_of_exact_solution(x_right_bound_, time_iter * tau_);

        local_result_.clear();

        ta_result = get_time_layer_result();

        free_part_.clear();

        local_result_.emplace_back(function_of_exact_solution(x_left_bound_, time_iter * tau_));
        local_result_.insert(local_result_.end(), ta_result.begin(), ta_result.end());
        local_result_.emplace_back(function_of_exact_solution(x_right_bound_, time_iter * tau_));

        results_.emplace_back(local_result_);
    }
    local_result_.clear();
}

std::vector<double> tridiagonal_matrix::get_time_layer_result() {
    auto n = free_part_.size();

    std::vector<double> result(n);

    std::vector<double> alpha(n - 1), beta(n - 1);

    double common_factor;

    alpha[n - 2] = -matrix_above_ / matrix_main_;
    beta[n - 2] = free_part_.back() / matrix_main_;

    for (auto iter = static_cast<long long int>(n) - 3; iter >= 0; iter--) {
        common_factor = 1.0 / (matrix_main_ + matrix_above_ * alpha[iter + 1]);
        alpha[iter] = -matrix_above_ * common_factor;
        beta[iter] = (free_part_[iter + 1] - beta[iter + 1] * matrix_above_) * common_factor;
    }

    result[0] = (free_part_[0] - matrix_above_ * beta[0]) / (matrix_main_ + matrix_above_ * alpha[0]);

    for (unsigned long iter = 1; iter < n; ++iter)
        result[iter] = alpha[iter - 1] * result[iter - 1] + beta[iter - 1];

    return result;
}

void tridiagonal_matrix::write_result() const {
    std::fstream result_file;
    result_file.open("../result/result.txt", std::ios::out | std::ios::trunc);
    result_file << "[[" << x_right_bound_ << "," << x_left_bound_ << ", " << h_num_ << "],"
                << "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],[";
    for (unsigned time_iter = 0; time_iter < time_layers_num_; ++time_iter) {
        result_file << "[";
        for (unsigned space_iter = 0; space_iter < h_num_; ++space_iter) {
            result_file << results_[time_iter][space_iter];
            if (space_iter + 1 != h_num_)
                result_file << ",";
        }
        result_file << "]";
        if (time_iter + 1 != time_layers_num_)
            result_file << ",";
    }
    result_file << "]]\n";
    result_file.close();
}

double tridiagonal_matrix::get_error() {
    double abs_exact_solution, abs_our_solution, results_difference;
    for (unsigned time_iter = 1; time_iter < time_layers_num_ - 1; ++time_iter)
        for (unsigned space_iter = 0; space_iter < h_num_; ++space_iter) {
            abs_exact_solution = std::abs(function_of_exact_solution(space_iter * h_, time_iter * tau_));
            abs_our_solution = std::abs(results_[time_iter][space_iter]);
            results_difference = std::abs(abs_exact_solution - abs_our_solution);
            if (results_difference > error_)
                error_ = results_difference;
        }
    return error_;
}

double tridiagonal_matrix::get_max_error() const {
    return tau_ + h_ * h_;
}
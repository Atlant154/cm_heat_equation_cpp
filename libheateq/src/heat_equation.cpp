#include "../include/heat_equation.h"

heat_equation::heat_equation(double (*heat_sources)(double, double),
                             double (*exact_solution)(double, double),
                             double diffusivity_coefficient,
                             uint32_t h_num,
                             uint32_t time_layers_num) :
      h_num_(h_num)
    , time_layers_num_(time_layers_num + 1)
    , a_(diffusivity_coefficient)
    , h_((x_right_bound_ - x_left_bound_) / h_num_)
    , tau_((time_right_bound_ - time_left_bound_) / time_layers_num_)
    , matrix_above_((-a_ * tau_) / pow(h_, 2))
    , matrix_main_(1.0 + 2.0 * a_ * tau_ / pow(h_, 2))
{
    std::vector<double> time_layer, ta_free_part, ta_result(h_num_ - 2);
    time_layer.reserve(time_layers_num_);
    ta_free_part.reserve(h_num_ - 2);
    ta_free_part.clear();

    results_.clear();
    results_.reserve(time_layers_num_);

    for (uint32_t iter = 0; iter < h_num_; ++iter)
        time_layer.emplace_back(exact_solution(iter * h_, time_left_bound_));

    results_.emplace_back(time_layer);

    double const courant_number = (a_ * tau_) / pow(h_, 2);

    for (uint32_t time_iter = 1; time_iter < time_layers_num_; ++time_iter)
    {
        for (uint32_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            ta_free_part.emplace_back(
                time_layer[space_iter]
                    + tau_ * heat_sources(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_));

        ta_free_part.front() += courant_number * exact_solution(x_left_bound_, time_left_bound_ + time_iter * tau_);
        ta_free_part.back() += courant_number * exact_solution(x_right_bound_, time_left_bound_ + time_iter * tau_);

        modified_thomas_alg(ta_free_part, ta_result);
        ta_free_part.clear();

        time_layer.clear();
        time_layer.emplace_back(exact_solution(x_left_bound_, time_iter * tau_));
        time_layer.insert(time_layer.end(), ta_result.begin(), ta_result.end());
        time_layer.emplace_back(exact_solution(x_right_bound_, time_iter * tau_));

        results_.emplace_back(time_layer);
    }
}

heat_equation::heat_equation(double (*heat_sources)(double, double),
                             double (*initial_time_layer)(double),
                             double (*left_bound)(double),
                             double (*right_bound)(double),
                             double diffusivity_coefficient,
                             uint32_t h_num,
                             uint32_t time_layers_num) :
      h_num_(h_num)
    , time_layers_num_(time_layers_num)
    , a_(diffusivity_coefficient)
    , h_((x_right_bound_ - x_left_bound_) / h_num_)
    , tau_((time_right_bound_ - time_left_bound_) / time_layers_num_)
    , matrix_above_((-a_ * tau_) / pow(h_, 2))
    , matrix_main_(1.0 + 2.0 * a_ * tau_ / pow(h_, 2))
{
    std::vector<double> time_layer, ta_free_part, ta_result(h_num_ - 2);
    time_layer.reserve(time_layers_num_);
    ta_free_part.reserve(h_num_ - 2);
    ta_free_part.clear();

    results_.clear();
    results_.reserve(time_layers_num_);

    for (std::size_t iter = 0; iter < h_num_; ++iter)
        time_layer.emplace_back(initial_time_layer(iter * h_));

    results_.emplace_back(time_layer);

    double const courant_number = (a_ * tau_) / pow(h_, 2);

    for (std::size_t time_iter = 1; time_iter < time_layers_num_; ++time_iter)
    {
        for (std::size_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            ta_free_part.emplace_back(
                time_layer[space_iter]
                    + tau_ * heat_sources(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_));

        ta_free_part.front() += courant_number * left_bound(time_left_bound_ + time_iter * tau_);
        ta_free_part.back() += courant_number * right_bound(time_left_bound_ + time_iter * tau_);

        modified_thomas_alg(ta_free_part, ta_result);
        ta_free_part.clear();

        time_layer.clear();
        time_layer.emplace_back(left_bound(time_left_bound_ + time_iter * tau_));
        std::copy(ta_result.begin(), ta_result.end(), std::back_inserter(time_layer));
        time_layer.emplace_back(right_bound(time_left_bound_ + time_iter * tau_));

        results_.emplace_back(time_layer);
    }
}

void heat_equation::modified_thomas_alg(std::vector<double> const & free_part, std::vector<double> & result)
{
    std::size_t n = result.size();
    std::vector<double> alpha(n - 1), beta(n - 1);

    double common_factor;

    alpha[n - 2] = -matrix_above_ / matrix_main_;
    beta[n - 2] = free_part.back() / matrix_main_;

    for (size_t iter = n - 3; iter > 0; iter--)
    {
        common_factor = 1.0 / (matrix_main_ + matrix_above_ * alpha[iter + 1]);
        alpha[iter] = -matrix_above_ * common_factor;
        beta[iter] = (free_part[iter + 1] - beta[iter + 1] * matrix_above_) * common_factor;
    }

    result[0] = (free_part[0] - matrix_above_ * beta[0]) / (matrix_main_ + matrix_above_ * alpha[0]);

    for (uint32_t iter = 1; iter < n; ++iter)
        result[iter] = alpha[iter - 1] * result[iter - 1] + beta[iter - 1];
}

void heat_equation::write_result(std::string const & path) const
{
    std::fstream result_file;
    result_file.open(path + "/tau:" + std::to_string(tau_) + "_h:" + std::to_string(h_) + "_result.txt",
                     std::ios::out | std::ios::trunc);
    result_file << "[[" << x_left_bound_ << "," << x_right_bound_ << ", " << h_num_ << "],"
                << "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],[";
    for (std::size_t time_iter = 0; time_iter < time_layers_num_; ++time_iter) {
        result_file << "[";
        for (std::size_t space_iter = 0; space_iter < h_num_; ++space_iter) {
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

void heat_equation::write_error_plot(double (*exact_solution)(double, double)) const
{
    double error;
    std::fstream result_file;
    result_file.open("../result/error.txt", std::ios::out | std::ios::trunc);
    result_file << "[[" << x_left_bound_ << "," << x_right_bound_ << ", " << h_num_ << "],"
                << "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],[";
    for (std::size_t time_iter = 0; time_iter < time_layers_num_; ++time_iter)
    {
        result_file << "[";
        for (std::size_t space_iter = 0; space_iter < h_num_; ++space_iter)
        {
            error = std::abs(exact_solution(space_iter * h_, time_iter * tau_) - results_[time_iter][space_iter]);
            result_file << error;
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

double heat_equation::get_error(double (*exact_solution)(double, double)) const
{
    double error = 0.0;
    double abs_appr, abs_exact;

    for (std::size_t time_iter = 1; time_iter < time_layers_num_; ++time_iter)
        for (std::size_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
        {
            abs_appr = std::abs(results_[time_iter][space_iter]);
            abs_exact = std::abs(exact_solution(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_));
            error += std::abs(abs_exact - abs_appr);
        }

    return error * tau_ * h_;
}
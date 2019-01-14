#include <heat_equation.h>

heat_equation::heat_equation(double_t (*heat_sources)(double_t, double_t),
                             double_t (*exact_solution)(double_t, double_t),
                             double_t diffusivity_coefficient,
                             uint32_t h_num,
                             uint32_t time_layers_num) :
      h_num_(h_num + 1)
    , time_layers_num_(time_layers_num + 1)
    , a_(diffusivity_coefficient)
    , h_((x_right_bound_ - x_left_bound_) / h_num_)
    , tau_((time_right_bound_ - time_left_bound_) / time_layers_num_)
    , matrix_above_((-a_ * tau_) / pow(h_, 2))
    , matrix_main_(1.0 + 2.0 * a_ * tau_ / pow(h_, 2))
{
    std::vector<double_t> time_layer, ta_free_part, ta_result(h_num_ - 2);
    time_layer.reserve(time_layers_num_);
    ta_free_part.reserve(h_num_ - 2);
    ta_free_part.clear();

    results_.clear();
    results_.reserve(time_layers_num_);

    for (uint32_t iter = 0; iter < h_num_; ++iter)
        time_layer.emplace_back(exact_solution(x_left_bound_ + iter * h_, time_left_bound_));

    results_.emplace_back(time_layer);

    double_t const courant_number = (a_ * tau_) / std::pow(h_, 2.);

    for (uint32_t time_iter = 1; time_iter < time_layers_num_; ++time_iter) {
        for (uint32_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            ta_free_part.emplace_back(
                    time_layer[space_iter]
                    + tau_ * heat_sources(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_));

        double_t const left_bound = exact_solution(x_left_bound_, time_left_bound_ + time_iter * tau_);
        double_t const right_bound = exact_solution(x_right_bound_, time_left_bound_ + time_iter * tau_);

        ta_free_part.front() += courant_number * left_bound;
        ta_free_part.back() += courant_number * right_bound;

        modified_thomas_alg(ta_free_part, ta_result);
        ta_free_part.clear();

        time_layer.clear();
        time_layer.emplace_back(left_bound);
        time_layer.insert(time_layer.end(), ta_result.begin(), ta_result.end());
        time_layer.emplace_back(right_bound);

        results_.emplace_back(time_layer);
    }
}

heat_equation::heat_equation(double_t (*heat_sources)(double_t, double_t),
                             double_t (*initial_time_layer)(double_t),
                             double_t (*left_bound)(double_t),
                             double_t (*right_bound)(double_t),
                             double_t diffusivity_coefficient,
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
    std::vector<double_t> time_layer, ta_free_part, ta_result(h_num_ - 2);
    time_layer.reserve(time_layers_num_);
    ta_free_part.reserve(h_num_ - 2);
    ta_free_part.clear();

    results_.clear();
    results_.reserve(time_layers_num_);

    for (uint32_t iter = 0; iter < h_num_; ++iter)
        time_layer.emplace_back(initial_time_layer(x_left_bound_ + iter * h_));

    results_.emplace_back(time_layer);

    double_t const courant_number = (a_ * tau_) / std::pow(h_, 2.);

    for (uint32_t time_iter = 1; time_iter < time_layers_num_; ++time_iter)
    {
        for (uint32_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            ta_free_part.emplace_back(
                time_layer[space_iter]
                    + tau_ * heat_sources(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_));

        double_t const left_bound_val = left_bound(time_left_bound_ + time_iter * tau_);
        double_t const right_bound_val = right_bound(time_left_bound_ + time_iter * tau_);

        ta_free_part.front() += courant_number * left_bound_val;
        ta_free_part.back() += courant_number * right_bound_val;

        modified_thomas_alg(ta_free_part, ta_result);
        ta_free_part.clear();

        time_layer.clear();
        time_layer.emplace_back(left_bound_val);
        time_layer.insert(time_layer.end(), ta_result.begin(), ta_result.end());
        time_layer.emplace_back(right_bound_val);

        results_.emplace_back(time_layer);
    }
}

void heat_equation::modified_thomas_alg(std::vector<double_t> const & free_part, std::vector<double_t> & result) {
    std::size_t n = result.size();
    std::vector<double_t> alpha(n - 1), beta(n - 1);

    double_t common_factor;

    alpha[n - 2] = -matrix_above_ / matrix_main_;
    beta[n - 2] = free_part.back() / matrix_main_;

    for (auto iter = n - 2; iter > 0; --iter) {
        common_factor = 1.0 / (matrix_main_ + matrix_above_ * alpha[iter]);
        alpha[iter - 1] = -matrix_above_ * common_factor;
        beta[iter - 1] = (free_part[iter] - beta[iter] * matrix_above_) * common_factor;
    }

    result[0] = (free_part[0] - matrix_above_ * beta[0]) / (matrix_main_ + matrix_above_ * alpha[0]);

    for (std::size_t iter = 1; iter < n; ++iter)
        result[iter] = alpha[iter - 1] * result[iter - 1] + beta[iter - 1];
}

void heat_equation::write_result(std::string const & path) const
{
    std::fstream result_file;
    result_file.open(path + "/result.txt",
                     std::ios::out | std::ios::trunc);
    result_file << "[[" << x_left_bound_ << "," << x_right_bound_ << ", " << h_num_ << "],"
                << "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],[";
    for (uint32_t time_iter = 0; time_iter < time_layers_num_; ++time_iter) {
        result_file << "[";
        for (uint32_t space_iter = 0; space_iter < h_num_; ++space_iter) {
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

void heat_equation::write_error_plot(double_t (*exact_solution)(double_t, double_t), std::string const & path) const
{
    double error;
    std::fstream result_file;
    result_file.open(path + "/error.txt", std::ios::out | std::ios::trunc);
    result_file << "[[" << x_left_bound_ << "," << x_right_bound_ << ", " << h_num_ << "],"
                << "[" << time_left_bound_ << "," << time_right_bound_ << ", " << time_layers_num_ << "],[";
    for (uint32_t time_iter = 0; time_iter < time_layers_num_; ++time_iter) {
        result_file << "[";
        for (uint32_t space_iter = 0; space_iter < h_num_; ++space_iter) {
            error = std::abs(exact_solution(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_) -
                             results_[time_iter][space_iter]);
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

double heat_equation::get_error(double_t (*exact_solution)(double_t, double_t)) const
{
    double_t error = 0.0;

    for (uint32_t time_iter = 1; time_iter < time_layers_num_; ++time_iter)
        for (uint32_t space_iter = 1; space_iter < h_num_ - 1; ++space_iter)
            error += std::abs(results_[time_iter][space_iter] -
                              exact_solution(x_left_bound_ + space_iter * h_, time_left_bound_ + time_iter * tau_)) *
                                      tau_ * h_;

    return error - std::numeric_limits<double_t>::epsilon() * h_num_ * time_layers_num_;
}

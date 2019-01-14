#include <heat_equation.h>
#include "gtest/gtest.h"

double_t const diffusivity_coefficient = 0.020;

inline double_t HeatSources(double_t const x, double_t const t) {
    return 2. * t - std::exp(x) + x - diffusivity_coefficient * ((-1.) * t * exp(x) - 12. * std::pow(x, 2.));
}

inline double_t ExactSolution(double_t const x, double_t const t) {
    return (-1.) * std::pow(x, 4.) + t * x + std::pow(t, 2.) - t * std::exp(x);
}

std::tuple<double_t, double_t> GetMethodError(uint32_t const first_h, uint32_t const first_tau, uint32_t factor) {
    uint32_t const second_h = first_h * factor;
    uint32_t const second_tau = first_tau * factor;

    heat_equation const first_split(HeatSources, ExactSolution, diffusivity_coefficient, first_h, first_tau);
    heat_equation const second_split(HeatSources, ExactSolution, diffusivity_coefficient, second_h, second_tau);

    double_t const first_error = first_split.get_error(ExactSolution);
    double_t const second_error = second_split.get_error(ExactSolution);

    return {first_error, second_error};
}

void HeatEquationErrorCorrection(double_t & second_error, uint32_t const factor) {
    second_error *= factor + std::pow(factor, 2);
}

TEST(Accuracy, LargeMesh) {
    uint32_t const first_h{10}, first_tau{10}, factor{2};
    auto[first_error, second_error] = GetMethodError(first_h, first_tau, factor);
    HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, second_error);
}

TEST(Accuracy, FineMesh) {
    uint32_t const first_h{10}, first_tau{10}, factor{500};
    auto[first_error, second_error] = GetMethodError(first_h, first_tau, factor);
    HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, second_error);
}

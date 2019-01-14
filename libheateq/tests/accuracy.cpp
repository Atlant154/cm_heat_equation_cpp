#include <heat_equation.h>
#include "gtest/gtest.h"

double const diffusivity_coefficient = 0.020;

inline double HeatSources(double_t const x, double_t const t) {
    return 2. * t - exp(x) + x - diffusivity_coefficient * ((-1.) * t * exp(x) - 12. * pow(x, 2));
}

inline double ExactSolution(double_t const x, double_t const t) {
    return (-1.) * pow(x, 4.) + t * x + pow(t, 2.) - t * exp(x);
}

std::tuple<double, double> GetMethodError(uint32_t const first_h, uint32_t const first_tau, uint32_t factor) {
    uint32_t const second_h = first_h * factor;
    uint32_t const second_tau = first_tau * factor;

    heat_equation first_split(HeatSources, ExactSolution, diffusivity_coefficient, first_h, first_tau);
    heat_equation second_split(HeatSources, ExactSolution, diffusivity_coefficient, second_h, second_tau);

    double_t first_error = first_split.get_error(ExactSolution);
    double_t second_error = second_split.get_error(ExactSolution);

    return {first_error, second_error};
}

void HeatEquationErrorCorrection(double_t & second_error, uint32_t factor) {
    second_error *= factor + std::pow(factor, 2);
}

TEST(Accuracy, LargeMesh) {
    uint32_t first_h{10}, first_tau{10}, factor{2};
    auto[first_error, second_error] = GetMethodError(first_h, first_tau, factor);
    HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, second_error);
}

TEST(Accuracy, FineMesh) {
    uint32_t first_h{10}, first_tau{10}, factor{500};
    auto[first_error, second_error] = GetMethodError(first_h, first_tau, factor);
    HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, second_error);
}
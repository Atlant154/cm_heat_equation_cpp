#include <heat_equation.h>
#include "gtest/gtest.h"

double_t const diffusivity_coefficient{0.020};

inline double_t HeatSources(double_t const x, double_t const t) {
    return 2. * t - std::exp(x) + x - diffusivity_coefficient * ((-1.) * t * std::exp(x) - 12. * std::pow(x, 2));
}

inline double_t ExactSolution(double_t const x, double_t const t) {
    return (-1.) * std::pow(x, 4) + t * x + std::pow(t, 2) - t * std::exp(x);
}

double_t HeatEquationApprError(uint32_t const h_num, uint32_t const tau_num){
    heat_equation const heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num);
    return heat_eq.get_error(ExactSolution);
}

double_t HeatEquationErrorCorrection(double_t const second_error, uint32_t const factor) {
    return second_error * factor + std::pow(factor, 2);
}

TEST(Accuracy, LargeMesh) {
    uint32_t const h_num{10}, tau_num{10}, factor{2};
    double_t const first_error{HeatEquationApprError(h_num, tau_num)};
    double_t const second_error{HeatEquationApprError(h_num * factor, tau_num * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

TEST(Accuracy, FineMesh) {
    uint32_t const h_num{10}, tau_num{10}, factor{500};
    double_t const first_error{HeatEquationApprError(h_num, tau_num)};
    double_t const second_error{HeatEquationApprError(h_num * factor, tau_num * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

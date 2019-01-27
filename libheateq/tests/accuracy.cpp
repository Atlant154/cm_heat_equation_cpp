#include <heat_equation.h>
#include "gtest/gtest.h"

double_t constexpr diffusivity_coefficient{0.010417};

inline double_t HeatSources(double_t const x, double_t const t) {
    return std::pow(x, 2.) * (std::pow(x, 2.)
                             - 12. * diffusivity_coefficient * t)
                             + std::exp(x) * (x * t * (diffusivity_coefficient * t - 2.)
                             + diffusivity_coefficient * std::pow(t, 2.) + 2. * t);
}

inline double_t ExactSolution(double_t const x, double_t const t) {
    return t * std::pow(x, 4.) - std::pow(t, 2.) * std::exp(x) * (x - 1.) + 1.;
}

inline double_t constexpr InitialTimeLayer(double_t const x) {
    return 1.;
}

inline double_t LeftBound(double_t const t) {
    return 1. + std::pow(t, 2.);
}

inline double_t RightBound(double_t const t) {
    return 1. + t;
}

double_t HeatEquationApprErrorFromExact(uint32_t const h_num, uint32_t const tau_num){
    heat_equation const heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num);
    return heat_eq.get_error(ExactSolution);
}

double_t HeatEquationApprErrorFromBounds(uint32_t const h_num, uint32_t const tau_num){
    heat_equation const heat_eq(HeatSources, InitialTimeLayer, LeftBound, RightBound, diffusivity_coefficient, h_num, tau_num);
    return heat_eq.get_error(ExactSolution);
}

double_t HeatEquationErrorCorrection(double_t const second_error, uint32_t const factor) {
    return second_error * std::pow(factor, 2.);
}

TEST(FromExactSolution, LargeMesh) {
    uint32_t const h_num{5}, tau_num{25}, factor{2};
    double_t const first_error{HeatEquationApprErrorFromExact(h_num, tau_num)};
    double_t const second_error{HeatEquationApprErrorFromExact(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

TEST(FromExactSolution, FineMesh) {
    uint32_t const h_num{5}, tau_num{25}, factor{100};
    double_t const first_error{HeatEquationApprErrorFromExact(h_num, tau_num)};
    double_t const second_error{HeatEquationApprErrorFromExact(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

TEST(FromBoundConditions, LargeMesh) {
    uint32_t const h_num{5}, tau_num{25}, factor{2};
    double_t const first_error{HeatEquationApprErrorFromBounds(h_num, tau_num)};
    double_t const second_error{HeatEquationApprErrorFromBounds(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

TEST(FromBoundConditions, FineMesh) {
    uint32_t const h_num{5}, tau_num{25}, factor{100};
    double_t const first_error{HeatEquationApprErrorFromBounds(h_num, tau_num)};
    double_t const second_error{HeatEquationApprErrorFromBounds(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error = HeatEquationErrorCorrection(second_error, factor);

    ASSERT_GE(first_error, corrected_error);
}

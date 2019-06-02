#include <HeatEquation.hpp>
#include "gtest/gtest.h"

double_t constexpr diffusivity_coefficient{0.010417};

std::function<double(double, double)> HeatSources = []( double X, double T ) -> double {
    return std::pow(X, 2.) * (std::pow(X, 2.) - 12. * diffusivity_coefficient * T)
           + std::exp(X) * (X * T * (diffusivity_coefficient * T - 2.) + diffusivity_coefficient * std::pow(T, 2.) + 2. * T);
};

std::function<double(double, double)> ExactSolution = []( double X, double T ) -> double {
    return T * std::pow( X, 4.0 ) - std::pow( T, 2.0 ) * std::exp( X ) * ( X - 1.0 ) + 1.0;
};

double_t HeatEquationApprErrorFromExact(uint32_t const h_num, uint32_t const tau_num){
    HeatEquation const heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num);
    return heat_eq.GetError(ExactSolution);
}

double_t HeatEquationErrorCorrection(double_t const second_error, uint32_t const factor) {
    return second_error * std::pow(factor, 2.);
}

TEST(FromExactSolution, LargeMesh) {
    uint32_t const h_num{5}, tau_num{h_num * h_num}, factor{2};
    double_t const first_error  {HeatEquationApprErrorFromExact(h_num, tau_num)};
    double_t const second_error {HeatEquationApprErrorFromExact(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error {HeatEquationErrorCorrection(second_error, factor)};

    ASSERT_GE(first_error, corrected_error);
}

TEST(FromExactSolution, FineMesh) {
    uint32_t const h_num{5}, tau_num{h_num * h_num}, factor{100};
    double_t const first_error  {HeatEquationApprErrorFromExact(h_num, tau_num)};
    double_t const second_error {HeatEquationApprErrorFromExact(h_num * factor, tau_num * factor * factor)};

    double_t corrected_error {HeatEquationErrorCorrection(second_error, factor)};

    ASSERT_GE(first_error, corrected_error);
}

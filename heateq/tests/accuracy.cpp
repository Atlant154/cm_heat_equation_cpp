#include <HeatEquation.hpp>

#include <gtest/gtest.h>

double constexpr diffusivity_coefficient{0.010417};

std::function<double(double, double)> HeatSources = []( double X, double T ) -> double {
    return std::pow(X, 2.) * (std::pow(X, 2.) - 12. * diffusivity_coefficient * T)
           + std::exp(X) * (X * T * (diffusivity_coefficient * T - 2.) + diffusivity_coefficient * std::pow(T, 2.) + 2. * T);
};

std::function<double(double, double)> ExactSolution = []( double X, double T ) -> double {
    return T * std::pow( X, 4.0 ) - std::pow( T, 2.0 ) * std::exp( X ) * ( X - 1.0 ) + 1.0;
};

double inline HeatEquationApprErrorFromExact( uint32_t const h_num, uint32_t const tau_num ){
    HeatEquation heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num, 7);
    return heat_eq.GetError(ExactSolution);
}

double inline HeatEquationErrorCorrection( double const second_error, uint32_t const factor) {
    return second_error * std::pow(factor, 2.);
}

TEST(FromExactSolution, LargeMesh) {
    uint32_t const h_num{ 5 }, tau_num{ h_num * h_num }, factor{ 2 };
    double const firstError  = HeatEquationApprErrorFromExact(h_num, tau_num);
    double const secondError = HeatEquationApprErrorFromExact(h_num * factor, tau_num * factor * factor);

    double const corrected_error = HeatEquationErrorCorrection(secondError, factor);

    EXPECT_NE( firstError, std::numeric_limits<decltype( firstError ) >::infinity() );
    EXPECT_NE( corrected_error, std::numeric_limits<decltype( corrected_error )>::infinity() );

    ASSERT_GE( firstError, corrected_error );
}

TEST(FromExactSolution, FineMesh) {
    uint32_t const h_num{ 5 }, tau_num{ h_num * h_num }, factor{ 100 };
    double const firstError  = HeatEquationApprErrorFromExact( h_num, tau_num );
    double const secondError = HeatEquationApprErrorFromExact( h_num * factor, tau_num * factor * factor );

    double const corrected_error = HeatEquationErrorCorrection(secondError, factor);

    EXPECT_NE( firstError, std::numeric_limits<decltype( firstError ) >::infinity() );
    EXPECT_NE( corrected_error, std::numeric_limits<decltype( corrected_error )>::infinity() );

    EXPECT_GE( firstError, corrected_error );
}

int32_t main( int32_t argc, char ** argv ) {
    testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
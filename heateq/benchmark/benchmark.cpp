#include <HeatEquation.hpp>
#include <benchmark/benchmark.h>

static void HeatEQ(benchmark::State &state) {
    double constexpr diffusivityCoefficient { 0.010417 };

    std::function<double(double, double)> HeatSources = []( double X, double T ) -> double
    {
        return std::pow(X, 2.) * (std::pow(X, 2.) - 12. * diffusivityCoefficient * T)
               + std::exp(X) * (X * T * (diffusivityCoefficient * T - 2.) + diffusivityCoefficient * std::pow(T, 2.) + 2. * T);
    };

    std::function<double(double, double)> ExactSolution = []( double X, double T ) -> double
    {
        return T * std::pow( X, 4.0 ) - std::pow( T, 2.0 ) * std::exp( X ) * ( X - 1.0 ) + 1.0;
    };

    for (auto _ : state)
        HeatEquation heq_eq(HeatSources, ExactSolution, diffusivityCoefficient, state.range(0), 100, 7);
}

BENCHMARK(HeatEQ)->RangeMultiplier(2)->Range(128, 1024<<10)->Iterations(20);
BENCHMARK_MAIN();

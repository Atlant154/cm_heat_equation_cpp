#include <heat_equation.h>
#include <benchmark/benchmark.h>

double const diffusivity_coefficient = 0.020;

inline double HeatSources(double x, double t) {
    return 2 * t - exp(x) + x - diffusivity_coefficient * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

inline double ExactSolution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void h_equation(benchmark::State & state)
{
    for (auto _ : state)
        heat_equation exmpl(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
}

BENCHMARK(h_equation)->RangeMultiplier(2)->Range(128, 1024<<7);

BENCHMARK_MAIN();
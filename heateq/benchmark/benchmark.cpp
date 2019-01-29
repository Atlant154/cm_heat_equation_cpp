#include <HeatEquation.hpp>
#include <benchmark/benchmark.h>

double_t constexpr diffusivity_coefficient{0.010417};

inline double_t HeatSources(double_t const x, double_t const t) noexcept {
    return std::pow(x, 2.) * (std::pow(x, 2.)
                             - 12. * diffusivity_coefficient * t)
                             + std::exp(x) * (x * t * (diffusivity_coefficient * t - 2.)
                             + diffusivity_coefficient * std::pow(t, 2.) + 2. * t);
}

inline double_t ExactSolution(double_t const x, double_t const t) noexcept {
    return t * std::pow(x, 4.) - std::pow(t, 2.) * std::exp(x) * (x - 1.) + 1.;
}

static void HeatEQ(benchmark::State &state) {
    for (auto _ : state)
        HeatEquation heq_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
}

BENCHMARK(HeatEQ)->RangeMultiplier(2)->Range(128, 1024<<7);

static void WriteJSON(benchmark::State & state){
    HeatEquation heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
    for(auto _ : state)
        heat_eq.WriteResultJSON();
}

BENCHMARK(WriteJSON)->Arg(256)->Arg(1024);

static void WriteOfstream(benchmark::State & state){
    HeatEquation heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
    for(auto _ : state)
        heat_eq.WriteResult();
}

BENCHMARK(WriteOfstream)->Arg(256)->Arg(1024);

BENCHMARK_MAIN();

#include <heat_equation.h>
#include <benchmark/benchmark.h>

double const diffusivity_coefficient{0.020};

inline double_t HeatSources(double_t const x, double_t const t) {
    return 2. * t - std::exp(x) + x - diffusivity_coefficient * ((-1.) * t * std::exp(x) - 12. * std::pow(x, 2));
}

inline double_t ExactSolution(double_t const x, double_t const t) {
    return (-1.) * std::pow(x, 4) + t * x + std::pow(t, 2) - t * std::exp(x);
}

static void HeatEQ(benchmark::State &state) {
    for (auto _ : state)
        heat_equation heq_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
}

BENCHMARK(HeatEQ)->RangeMultiplier(2)->Range(128, 1024<<7);

static void WriteJSON(benchmark::State & state){
    heat_equation heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
    for(auto _ : state){
        heat_eq.write_result_json();
    }
}

BENCHMARK(WriteJSON)->Arg(256)->Arg(1024);

static void Write(benchmark::State & state){
    heat_equation heat_eq(HeatSources, ExactSolution, diffusivity_coefficient, state.range(0), 100);
    for(auto _ : state){
        heat_eq.write_result();
    }
}

BENCHMARK(Write)->Arg(256)->Arg(1024);

BENCHMARK_MAIN();
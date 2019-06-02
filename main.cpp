#include <iostream>

#include "CLI/CLI11.hpp"
#include <HeatEquation.hpp>

double_t constexpr diffusivity_coefficient{0.010417};

inline double_t ExactSolution(double_t const x, double_t const t) noexcept {
    return t * std::pow(x, 4.) - std::pow(t, 2.) * std::exp(x) * (x - 1.) + 1.;
}

int32_t main(int32_t argc, char * argv[]) {
    CLI::App application{"Application is an example of solving the heat equation"};

    uint32_t tau_num{0};
    application.add_option("-t,--tau", tau_num, "The number of time layers")->required();

    uint32_t h_num{0};
    application.add_option("-s,--splits", h_num, "The number of of spatial splits")->required();

    bool write_to_file{false};
    application.add_flag("-w,--write", write_to_file, "Write results to file");

    std::string output_dir{"."};
    application.add_option("-o,--output-path", output_dir,
                           "The output directory path, i.e. ~/projects/visualization. Default: ./")->required(false);

    CLI11_PARSE(application, argc, argv);

    std::function<double(double, double)> HeatSources = []( double X, double T ) -> double {
        return std::pow(X, 2.) * (std::pow(X, 2.) - 12. * diffusivity_coefficient * T)
            + std::exp(X) * (X * T * (diffusivity_coefficient * T - 2.) + diffusivity_coefficient * std::pow(T, 2.) + 2. * T);
    };

    std::function<double(double, double)> ExactSolution = []( double X, double T ) -> double {
        return T * std::pow( X, 4.0 ) - std::pow( T, 2.0 ) * std::exp( X ) * ( X - 1.0 ) + 1.0;
    };

    HeatEquation hq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num);

    double_t const error = hq.GetError(ExactSolution);
    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;

//    if (write_to_file) {
//        hq.WriteResultJSON(output_dir);
//        hq.WriteErrorJSON(ExactSolution, output_dir);
//        hq.WriteExactSolutionJSON(ExactSolution, output_dir);
//    }

    return EXIT_SUCCESS;
}
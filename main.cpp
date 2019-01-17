#include <iostream>

#include <CLI/CLI.hpp>
#include <heat_equation.h>

double_t const diffusivity_coefficient{0.020};

inline double_t HeatSources(double_t const x, double_t const t) {
    return 2. * t - std::exp(x) + x - diffusivity_coefficient * ((-1.) * t * std::exp(x) - 12. * std::pow(x, 2));
}

inline double_t ExactSolution(double_t const x, double_t const t) {
    return (-1.) * std::pow(x, 4) + t * x + std::pow(t, 2) - t * std::exp(x);
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

    heat_equation hq(HeatSources, ExactSolution, diffusivity_coefficient, h_num, tau_num);

    double_t const error = hq.get_error(ExactSolution);
    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;

    if (write_to_file) {
        hq.write_result(output_dir);
        hq.write_error_plot(ExactSolution, output_dir);
    }

    return EXIT_SUCCESS;
}
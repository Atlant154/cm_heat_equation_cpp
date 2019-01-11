#include <iostream>
#include <getopt.h>

#include <heat_equation.h>

double const diffusivity_coefficient = 0.020;
std::string g_path = ".";
bool g_write = false;
bool g_help = false;
uint32_t g_tau_num = 0;
uint32_t g_h_num = 0;

double function_of_heat_sources(double x, double t) {
    return 2 * t - exp(x) + x - diffusivity_coefficient * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double function_of_exact_solution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void show_usage()
{
    std::cout << "Heat equation tool keys usage:"
        << '\n' << "-t --tau         [number]      [required] The number of time layers."
        << '\n' << "-s --splits      [number]      [required] The number of of spatial splits."
        << '\n' << "-w --write                     [optional] Write results to file."
        << '\n' << "-o --output-path [path]        [optional] The output directory path, i.e. ~/projects/result. Default: ./."
        << '\n' << "-h --help                      [optional] Show this message."
        << '\n';
}

bool parse_arguments(int argc, char ** argv) {
    static struct option long_options[] =
            {{"tau",         required_argument, nullptr, 't'},
             {"splits",      required_argument, nullptr, 's'},
             {"output-path", required_argument, nullptr, 'o'},
             {"write",       no_argument,       nullptr, 'w'},
             {"help",        no_argument,       nullptr, 'h'}
            };

    int32_t opt;
    while ((opt = getopt(argc, argv, "t:s:o:h:w")) != -1) {
        switch (opt) {
            case 't': {
                int64_t const tau_candidate = std::stoi(optarg);
                g_tau_num = tau_candidate > 0 ? (uint32_t) tau_candidate : 0;
            }
                break;
            case 's': {
                int64_t const h_candidate = std::stoi(optarg);
                g_h_num = h_candidate > 0 ? (uint32_t) h_candidate : 0;
            }
                break;
            case 'w':
                g_write = true;
                break;
            case 'h':
                g_help = true;
                break;
            default:
                g_help = true;
        }
    }
    return g_h_num != 0 && g_tau_num != 0 && !g_help;
}

int main(int argc, char ** argv) {
    if (!parse_arguments(argc, argv)) {
        show_usage();
        return EXIT_FAILURE;
    }

    heat_equation
            test_class_member(function_of_heat_sources, function_of_exact_solution, diffusivity_coefficient, g_h_num,
                              g_tau_num);

    double error = test_class_member.get_error(function_of_exact_solution);
    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;

    if (g_write) {
        test_class_member.write_result(g_path);
        test_class_member.write_error_plot(function_of_exact_solution);
    }

    return EXIT_SUCCESS;
}
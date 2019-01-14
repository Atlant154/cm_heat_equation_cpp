#include <iostream>
#include <getopt.h>

#include <heat_equation.h>

double const diffusivity_coefficient = 0.020;

std::string g_outputPath = ".";
bool g_write = false;
bool g_help = false;
uint32_t g_tauNum = 0;
uint32_t g_hNum = 0;

inline double HeatSources(double x, double t) {
    return 2 * t - exp(x) + x - diffusivity_coefficient * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

inline double ExactSolution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

void show_usage()
{
    std::cout << "Heat equation tool keys usage:"
        << '\n' << "-t --tau         [number]      [required] The number of time layers."
        << '\n' << "-s --splits      [number]      [required] The number of of spatial splits."
        << '\n' << "-w --write                     [optional] Write results to file."
        << '\n' << "-o --output-path [path]        [optional] The output directory path, i.e. ~/projects/visualization. Default: ./."
        << '\n' << "-h --help                      [optional] Show this message."
        << '\n';
}

bool parse_arguments(int argc, char ** argv) {
    static struct option long_options[] =
            {{"tau",         required_argument, nullptr, 't'},
             {"splits",      required_argument, nullptr, 's'},
             {"output-path", required_argument, nullptr, 'o'},
             {"write",       no_argument,       nullptr, 'w'},
             {"help",        no_argument,       nullptr, 'h'},
             {nullptr,       0,                 nullptr, 0  }
            };

    int32_t opt;
    int32_t opt_index;
    while ((opt = getopt_long(argc, argv, "t:s:o:h:w", long_options, &opt_index)) != -1) {
        switch (opt) {
            case 't': {
                int64_t const tau_candidate = std::stoi(optarg);
                g_tauNum = tau_candidate > 0 ? (uint32_t) tau_candidate : 0;
            }
                break;
            case 's': {
                int64_t const h_candidate = std::stoi(optarg);
                g_hNum = h_candidate > 0 ? (uint32_t) h_candidate : 0;
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
    return g_hNum != 0 && g_tauNum != 0 && !g_help;
}

int main(int argc, char ** argv) {
    if (!parse_arguments(argc, argv)) {
        show_usage();
        return EXIT_FAILURE;
    }

    heat_equation
            test_class_member(HeatSources, ExactSolution, diffusivity_coefficient, g_hNum,
                              g_tauNum);

    double error = test_class_member.get_error(ExactSolution);
    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;
    std::cout << "Theoretical error: " << 1. / (g_hNum * g_hNum + g_tauNum) << "." << std::endl;

    if (g_write) {
        test_class_member.write_result(g_outputPath);
        test_class_member.write_error_plot(ExactSolution, g_outputPath);
    }

    return EXIT_SUCCESS;
}
#include <iostream>

#include "libheateq/include/heat_equation.h"

double const diffusivity_coefficient = 0.020;

double function_of_heat_sources(double x, double t) {
    return 2 * t - exp(x) + x - diffusivity_coefficient * ((-1) * t * exp(x) - 12 * pow(x, 2));
}

double function_of_exact_solution(double x, double t) {
    return (-1) * pow(x, 4) + t * x + pow(t, 2) - t * exp(x);
}

int main() {
    uint32_t h_num = 0;
    while (true) {
        std::cout << "Hi. Enter here the number of partitions of the object: ";
        std::cin >> h_num;
        if (h_num > 0)
            break;
    }

    uint32_t time_layers_num = 0;
    while (true) {
        std::cout << "Put number of time layers here: ";
        std::cin >> time_layers_num;
        if (time_layers_num > 0)
            break;
    }

    heat_equation
        test_class_member(function_of_heat_sources, function_of_exact_solution, diffusivity_coefficient, h_num, time_layers_num);

    double error = test_class_member.get_error(function_of_exact_solution);
    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;

    test_class_member.write_result();
    test_class_member.write_error_plot(function_of_exact_solution);

    return EXIT_SUCCESS;
}
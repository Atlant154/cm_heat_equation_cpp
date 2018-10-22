#include <iostream>

#include "include/heat_equation.h"

int main() {
    unsigned int h_num = 0;

    while (true) {
        std::cout << "Hi. Enter here the number of partitions of the object:" << std::endl;
        std::cin >> h_num;
        if (h_num > 0)
            break;
    }
    unsigned int time_layers_num = 0;

    while (true) {
        std::cout << "Put number of time layers here:" << std::endl;
        std::cin >> time_layers_num;
        if (time_layers_num > 0)
            break;
    }

    heat_equation test_class_member(h_num, time_layers_num);
    test_class_member.get_result();
    test_class_member.write_result();

    double error = test_class_member.get_error();
    double theoretical_error = test_class_member.get_max_error();

    std::cout << "Error resulting from the calculation: " << error << "." << std::endl;
    std::cout << "Theoretical error: " << theoretical_error << ". " << std::endl;

    if (error < theoretical_error)
        std::cout << "[Success] Real error less than theoretical!" << std::endl;
    else
        std::cout << "[Warning] Real error greater than theoretical!" << std::endl;

    return 0;
}
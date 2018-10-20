#include <iostream>

#include "include/tridiagonal_matrix.h"

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
    tridiagonal_matrix test_class_member(h_num, time_layers_num);
    test_class_member.get_result();
    test_class_member.write_result();
    std::cout << "Error is: " << test_class_member.get_error() << "." << std::endl;
    std::cout << "Expected error: " << test_class_member.get_max_error() << ". " << std::endl;
    return 0;
}
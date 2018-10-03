#include <iostream>

#include "../include/tridiagonal_matrix.h"

int main(){
    unsigned int h_num = 0;
    while(true){
        std::cout << "Hi. Enter here the number of partitions of the object:" << std::endl;
        std::cin >> h_num;
        if(h_num > 0)
            break;
    }
    unsigned int time_layers_num = 0;
    while(true){
        std::cout << "Put number of time layers here:" << std::endl;
        std::cin >> time_layers_num;
        if(time_layers_num > 0)
            break;
    }
    tridiagonal_matrix test_class_member(h_num, time_layers_num);
    test_class_member.get_result_();
    test_class_member.write_result();
    double error = test_class_member.get_error();
    double max_eror = test_class_member.get_max_error();
    std::cout << "Error is " << error << ".\n";
    std::cout << "Maximum treoretical error is " << max_eror << "." << std::endl;
    if(error > max_eror)
        std::cout << "(WARNING) Error greater than theoretical." << std::endl;
    else
        std::cout << "Error less than theretical." << std::endl;
    return 0;
}
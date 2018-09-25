#include <iostream>

#include "../include/tridiagonal_matrix.h"

int main(){
    long unsigned int h_num = 0;
    while(true){
        std::cout << "Hi. Enter here the number of partitions of the object:" << std::endl;
        std::cin >> h_num;
        if(h_num > 0)
            break
    }
    long unsigned int time_layers_num = 0;
    while(true){
        std::cout << "Put number of time layers here:" << std::endl;
        std::cin >> time_layers_num;
        if(time_layers_num > 0)
            break
    }
    return 0;
}
#include <iostream>
#include <fstream>

#include "../include/tridiagonal_matrix.h"

#define RESULT_PATH "../result/result.txt"

void write_to_file(std::vector<double> results, long unsigned int h_num, long unsigned int time_layers_num);

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
    //Test data:
    tridiagonal_matrix test_class_member(h_num, time_layers_num);
    std::vector<double> results = test_class_member.get_result_();
    write_to_file(results, h_num, time_layers_num);
    /* TODO make calculation great again */
    return 0;
}

void write_to_file(std::vector<double> results, long unsigned int h_num, long unsigned int time_layers_num)
{
    std::ofstream result_file;
    result_file.open(RESULT_PATH, std::ios::out | std::ios::trunc);
    /*TODO fix the getting bounds: */
    result_file << "[[" << 0 << ", " << 1 << ", " << (double)1/(h_num) << "],"
    << "[" << 0 << ", " << 1 << ", " << (double)1/(time_layers_num) << "],[";
    for(long unsigned int time_iter = 0; time_iter < time_layers_num; ++time_iter)
    {
        result_file << "[";
        for(long unsigned int h_iter = 0; h_iter < h_num; ++h_iter) {
            result_file << results[(time_iter) * (time_layers_num) + h_iter];
            std::cout << "H iter: " << h_iter << ". H num: " << h_num << "." << std::endl;
            if(h_iter + 1 != h_num)
                result_file << ",";
        }
        result_file << "]";
        if(time_iter + 1 != time_layers_num)
            result_file << ",";
    }
    result_file << "]]";
    result_file.close();
}
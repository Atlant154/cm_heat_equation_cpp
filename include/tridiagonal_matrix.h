#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#define X_LEFT_BOUND 0.0
#define X_RIGHT_BOUND 1.0
#define TIME_LEFT_BOUND 0.0
#define TIME_RIGHT_BOUND 1.0

#include<vector>
class tridiagonal_matrix {
public:	
	tridiagonal_matrix() {};
	~tridiagonal_matrix() {};
	
	tridiagonal_matrix(const unsigned int h_num, const unsigned int time_layers_num);
	
	long double get_h();
	
	std::vector<long double> TDMA();
private:
	void set_h(const unsigned int h_num);
	
	void set_tau(const unsigned int time_layers_num);
	
	void set_above();
	
	void set_main();
	
	void set_below();
	
	long double h = 0;
	unsigned int n = 0;
	long double tau = 0;
	std::vector<long double> above = {0};
	std::vector<long double> main = {0};
	std::vector<long double> below = {0};
	std::vector<long double> free = {0};
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H
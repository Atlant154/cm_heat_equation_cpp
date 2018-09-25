#ifndef PROJECT_TRIDIAGONAL_MATRIX_H
#define PROJECT_TRIDIAGONAL_MATRIX_H

#define X_LEFT_BOUND 0
#define X_RIGHT_BOUND 1
#define TIME_LEFT_BOUND 0
#define TIME_RIGHT_BOUND 1

#include<vector>
class tridiagonal_matrix {
public:	
	tridiagonal_matrix() {}
	~tridiagonal_matrix() {}
	
	tridiagonal_matrix(const long long h_num, const long long time_layers_num) {
		set_h(h_num);
		set_tau(time_layers_num);
		n = h_num;
		set_above();
		set_main();
		set_below()
	}
	
	long long get_h() {
		return h;
	};
private:
	void set_h(const long long h_num) {
		h = (X_RIGHT_BOUND - X_LEFT_BOUND) / h_num;
	}
	
	void set_tau(const long long time_layers_num) {
		tau = (TIME_RIGHT_BOUND - TIME_LEFT_BOUND) / time_layers_num;
	}
	
		void set_above() {
		above.reserve(n-1);
		for(int i=0; i < n-1; ++i) {
			above[i] = i;
		}
	}
	
	void set_main() {
		main.reserve(n)
		for(int i=0; i < n; ++i) {
			main[i] = 2*i*i;
		}		
	}
	
	void set_below() {
		below.reserve(n-1)
		for(int i=0; i < n-1; ++i) {
			below[i] = i;
		}		
	}
	
	long double h = 0;
	long long n = 0;
	long double tau = 0;
	std::vector<long double> above = {0};
	std::vector<long double> main = {0};
	std::vector<long double> below = {0};
	std::vector<long double> free = {0};
};


#endif //PROJECT_TRIDIAGONAL_MATRIX_H
#include "../include/tridiagonal_matrix.h"

tridiagonal_matrix::tridiagonal_matrix(const unsigned int h_num, const unsigned int time_layers_num) {
	set_h(h_num);
	set_tau(time_layers_num);
	n = h_num;
	set_above();
	set_main();
	set_below();
};

long double tridiagonal_matrix::get_h() {
	return h;
};

void tridiagonal_matrix::set_h(const unsigned int h_num) {
	h = (X_RIGHT_BOUND - X_LEFT_BOUND) / h_num;
};

void tridiagonal_matrix::set_tau(const unsigned int time_layers_num) {
	tau = (TIME_RIGHT_BOUND - TIME_LEFT_BOUND) / time_layers_num;
};

void tridiagonal_matrix::set_above() {
	above.reserve(n-1);
	for(int i=0; i < n-1; ++i) {
		above[i] = i;
	}
};

void tridiagonal_matrix::set_main() {
	main.reserve(n);
		for(int i=0; i < n; ++i) {
		main[i] = 2*i*i;
	}		
};
	
void tridiagonal_matrix::set_below() {
	below.reserve(n-1);
	for(int i=0; i < n-1; ++i) {
		below[i] = i;
	}		
};

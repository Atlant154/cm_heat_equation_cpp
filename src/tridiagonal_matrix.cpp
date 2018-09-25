#include "../include/tridiagonal_matrix.h"

tridiagonal_matrix::tridiagonal_matrix(const unsigned int h_num, const unsigned int time_layers_num) {
	set_h(h_num);
	set_tau(time_layers_num);
	n = h_num;
	set_above();
	set_main();
	set_below();
	
	set_a(0.022);
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

void set_a(long double alpha) {
	a = alpha;
};

void tridiagonal_matrix::set_above() {
	above.reserve(n-1);
};

void tridiagonal_matrix::set_main() {
	main.reserve(n);
	}		
};
	
void tridiagonal_matrix::set_below() {
	below.reserve(n-1);	
};

long double f(long double x, long double t) {
	return 1;
};

long double u(long double x, long double t) {
	return 1;
};

long double mu_0(long double x) {
	return 1;
}

long double mu_1(long double t) {
	return 1;
};
	
long double mu_2(long double t) {
	return 1;
};

std::vector<long double> tridiagonal_matrix::TDMA() {
	std::vector<long double> w;
	w.reserve(n-1);
	std::vector<long double> g;
	g.reserve(n);
	std::vector<long double> p;
	p.reserve(n);
	w[0] = below[0]/main[0];
	g[0] = free[0]/main[0];
	
	for(unsigned int i = 1; i < n-1; ++i) {
		w[i] = below[i]/(main[i] - above[i-1]*w[i-1]);
		g[i] = (free[i] - above[i-1]*g[i-1])/(main[i] - above[i-1]*w[i-1]);
	}
	g[n-1] = (free[n-1] - above[n-2]*g[n-2])/(main[n-1] - above[n-2]*w[n-2]);
	
	p[n-1] = g[n-1];
	for(unsigned int j = n-1; j > 0; --j) {
		p[j-1] = g[j-1] - w[j-1]*p[j];
	}
	return p;
};
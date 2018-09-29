#include "../include/tridiagonal_matrix.h"

tridiagonal_matrix::tridiagonal_matrix(unsigned int h_num, unsigned int time_layers_num) : h_num(h_num), time_layers_num(time_layers_num)
{
	set_h();
	set_tau();
	n = h_num;
	set_above();
	set_main();
	set_below();
}

long double tridiagonal_matrix::get_h_() {
	return h;
}

long double tridiagonal_matrix::get_tau_() {
    return tau;
}

void tridiagonal_matrix::set_h() {
	h = (x_right_bound - x_right_bound) / h_num;
}

void tridiagonal_matrix::set_tau() {
	tau = (time_right_bound - time_left_bound) / time_layers_num;
}

void tridiagonal_matrix::set_above() {
	above.reserve(n-1);
}

void tridiagonal_matrix::set_main() {
	main.reserve(n);
}

void tridiagonal_matrix::set_below() {
	below.reserve(n-1);
}

long double tridiagonal_matrix::f(long double x, long double t) {
	return 1;
}

long double tridiagonal_matrix::u(long double x, long double t) {
	return 1;
}

long double tridiagonal_matrix::mu_0(long double x) {
	return 1;
}

long double tridiagonal_matrix::mu_1(long double t) {
	return 1;
}

long double tridiagonal_matrix::mu_2(long double t) {
	return 1;
}

void tridiagonal_matrix::set_u() {
	/*TODO correct the logic*/
	u.reserve((h_num + 1) * time_layers_num)
	for(unsigned int j = 0; j < h_num + 1; ++j) {
		u[j] = mu_1(tau*j)
		u[n*h_num + j] = mu_2[tau*j]
	}

	for(unsigned int i = 0; i < time_layers_num; ++i) {
		u[h_num*i] = mu_0(h*) // WIP
	}
}

std::vector<long double> tridiagonal_matrix::get_result_()
{
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
}

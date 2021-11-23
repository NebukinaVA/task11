#include <iostream>
#include "Pendulum.h"

int main()
{
	/*
	    x0 = _x0;
		v0 = _v0;
		u0 = _u0;
		g = _g;
		L = _L;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
*/
	int i;
	Pendulum equation(0.0, 1.0, 1.0, 9.8, 5.0, 0.0001, 1000, 1e-6, 0.7, 1e-4);
	equation.calculate_w_error();
	//	equation.calculate();
	std::cout << "RK3 METHOD\n" << equation << std::endl;
	std::cin >> i;
	return 0;
}
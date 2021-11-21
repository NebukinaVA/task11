#pragma once
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <iterator>
#include <utility>
//Jj8kxf
// u'' + g/L*sin(u)

class Pendulum
{
private:
	long double u0, v0, x0;
	long double g, L;
	long double h, eps, xmax, prec;
	long int n;                   // number of steps
	std::vector<long double> arg; //x
//	std::vector<long double> ures; //v
//	std::vector<long double> vres; //u
	std::vector<std::pair<long double, long double>> res; // (u, v)
	std::vector<std::pair<long double, long double>> reswcap; //(u, v) with cap
	std::vector<long double> steps; //h
	std::vector<long double> ss;    //S*
//??? is it needed ?
//	std::vector<long double> exres; // u - exact result
	std::vector<long int> hinc;  // total step increases
	std::vector<long int> hdec;  // total step decreases
	long double func1(long double v)
	{
		return v;
	}
	long double func2(long double u)
	{
		return (-g * sin(u) / L);
	}
public:
	Pendulum(long double _x0, long double _u0, long double _v0, long double _g, long double _L, long double _h, long int _n, long double _eps, long double _xmax, long double _prec)
	{
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
	}
	std::pair<long double, long double> RK4(long double xn, long double un, long double vn, long double h)
	{
		// first u'=func1(v)
		long double k1 = func1(vn);
		long double k2 = func1(vn + h * k1 / 4.0);
		long double k3 = func1(vn + h * k2 / 2.0);
		long double k4 = func1(vn + h * (k1 - 2.0 * k2 + 2.0 * k3));
		un += h * (k1 + 4.0 * k3 + k4) / 6.0;
		// second v'= func2(u)
		k1 = func2(un);
		k2 = func2(un + h * k1 / 4.0);
		k3 = func2(un + h * k2 / 2.0);
		k4 = func2(un + h * (k1 - 2.0 * k2 + 2.0 * k3));
		vn += h * (k1 + 4.0 * k3 + k4) / 6.0;
		return std::make_pair(un, vn);
	}
	std::vector<std::pair<long double, long double>> calculate()
	{
	//	exres.push_back(v0);
		arg.push_back(x0);
		auto result = std::make_pair(u0, v0);
		res.push_back(result);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(h);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		long double xn = x0;
		long double un = u0;
		long double vn = v0;
		long int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn < xmax))
			{
				break;
			}
			else {
				if ((xn + h) > xmax)
				{
					while (((xn + h) > xmax) && (xn < (xmax - prec)))
					{
						h /= 2.0;
					}
					result = RK4(xn, result.first, result.second, h);
					xn += h;
				}
				else
				{
					result = RK4(xn, result.first, result.second, h);
					xn += h;
				}
				arg.insert(arg.begin() + i + 1, xn);
				res.insert(res.begin() + i + 1, std::make_pair(un, vn));
 				steps.insert(steps.begin() + i + 1, h);
			//	exres.insert(exres.begin() + i + 1, ExactSolution(xn));
				hinc.insert(hinc.begin() + i + 1, 0);
				hdec.insert(hdec.begin() + i + 1, 0);
				reswcap.insert(reswcap.begin() + i + 1, std::make_pair(0.0, 0.0));
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		return res;
	}
	std::vector<std::pair<long double, long double>> calculate_w_error()
	{
	//	exres.push_back(v0);
		auto result = std::make_pair(u0, v0);
		ss.push_back(0.0);
		arg.push_back(x0);
		res.push_back(result);
		steps.push_back(0.0);
		steps.push_back(h);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		hinc.push_back(0);
		hdec.push_back(0);
		hinc.push_back(0);
		hdec.push_back(0);
		long double xn = x0;
		long double vn = v0;
		long double xhalf = x0;
//		long double uhalf, vhalf, uwithcap, vwithcap, unext, vnext;
		auto half = result;
		auto cap = result;
		auto next = result;
		long double S, S1, S2;
		long int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn < xmax))
			{
				break;
			}
			else
			{
				
				half = RK4(xn, result.first, result.second, h / 2.0);
				xhalf = xn + h / 2.0;
			    cap = RK4(xhalf, half.first, half.second, h / 2.0);
				next = RK4(xn, result.first, result.second, h);
				S1 = (cap.first - next.first) / 15.0;
				S2 = (cap.second - next.second) / 15.0;
				S = sqrt(pow(S1, 2) + pow(S2, 2)); // норма погрешности
				if ((S >= (eps / 32.0)) && (S <= eps))
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
					//	vn = RK4(xn, vn, h);
						result = RK4(xn, result.first, result.second, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i + 1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i + 1]));
					}
					else
					{
						result = next;
						xn += h;
						hinc.insert(hinc.begin() + i + 2, hinc[i + 1]);
						hdec.insert(hdec.begin() + i + 2, hdec[i + 1]);
					}

				}
				else if (S < (eps / 16.0))
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i + 1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i + 1]));
					}
					else {
						result = next;
						xn += h;
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 2, ++(hinc[i + 1]));
						hdec.insert(hdec.begin() + i + 2, hdec[i + 1]);
					}
				}
				else if (S > eps)
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i + 1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i + 1]));
					}
					else {
						h = h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						xn += h;
						hinc.insert(hinc.begin() + i + 2, hinc[i + 1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i + 1]));
					}
				}
				i++;
				S *= 8.0;
				reswcap.insert(reswcap.begin() + i, cap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				res.insert(res.begin() + i, result);
				steps.insert(steps.begin() + i + 1, h);
			//	exres.insert(exres.begin() + i, ExactSolution(xn));
			}
		}
		return res;
	}
	long double ExactSolution(long double x)
	{
	//	return (v0*exp(2.0*x));
		return 0;
	}
	friend std::ostream & operator<<(std::ostream &out, Pendulum &vc)
	{
		if (vc.res.empty())
			out << "There are no calculated results yet.";
		else {
			out << "n  " << " h n-1  " << "    x    " << "       S*      " << "inc " << "dec" << std::endl;
		//	long double globerr;
			for (long int i = 0; i < vc.res.size(); i++)
			{
			//	globerr = abs(vc.exres[i] - vc.res[i]);
				out << i << "  " << vc.steps[i] << "      " << vc.arg[i] << "    " << vc.ss[i] << "    " << vc.hinc[i] << "  " << vc.hdec[i] << std::endl;
			}
		}
		return out;
	}

};
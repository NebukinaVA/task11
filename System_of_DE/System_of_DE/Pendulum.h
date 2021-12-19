#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

// u'' + g/L*sin(u)

class Pendulum
{
private:
	double u0, v0, x0;
	double g, L;
	double h, eps, xmax, prec;
	int n, N;                   // number of steps, total steps
	std::vector<double> arg; //x
	std::vector<double> ures; //v
	std::vector<double> vres; //u
	double Xn, Un, Vn;
	int inc = 0;
	int dec = 0; // total inc dec
	int Smin = 0, Smax = 0; // number of string with Smin, Smax
	int hmax, hmin; //number of string with hmax, hmin
//	std::vector<std::pair<double, double>> res; // (u, v)
	std::vector<std::pair<double, double>> reswcap; //(u, v) with cap
	std::vector<double> steps; //h
	std::vector<double> ss;    //S*
	std::vector<int> hinc;  // total step increases
	std::vector<int> hdec;  // total step decreases
	double func1(double x, double u, double v) // u' = v
	{
		return v;
	}
	double func2(double x, double u, double v) // v'= u'' = -g*sinu / L
	{
		return (-g * sin(u) / L);
	}
public:
	Pendulum(double _x0, double _u0, double _v0, double _g, double _L, double _h, int _n, double _eps, double _xmax, double _prec)
	{
		x0 = _x0;
		u0 = _u0;
		v0 = _v0;
		g = _g;
		L = _L;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
	}
	std::pair<double, double> RK4(double xn, double un, double vn, double h)
	{
		/*
        // first u'=func1(v)
		double k1 = func1(vn);
		double k2 = func1(vn + h * k1 / 4.0);
		double k3 = func1(vn + h * k2 / 2.0);
		double k4 = func1(vn + h * (k1 - 2.0 * k2 + 2.0 * k3));
		un += h * (k1 + 4.0 * k3 + k4) / 6.0;
		// second v'= func2(u)
		k1 = func2(un);
		k2 = func2(un + h * k1 / 4.0);
		k3 = func2(un + h * k2 / 2.0);
		k4 = func2(un + h * (k1 - 2.0 * k2 + 2.0 * k3));
		vn += h * (k1 + 4.0 * k3 + k4) / 6.0;
		*/
		// u' = v
		// v'= -g*sinu / L
		double ku1, ku2, ku3, ku4;  // for un
		double kv1, kv2, kv3, kv4;  // for vn

		ku1 = func1(xn, un, vn);
		kv1 = func2(xn, un, vn);

		ku2 = func1(xn + h / 4.0, un + (h / 4.0)*ku1, vn + (h / 4.0)*kv1);
		kv2 = func2(xn + h / 4.0, un + (h / 4.0)*ku1, vn + (h / 4.0)*kv1);
		
		ku3 = func1(xn + h / 2.0, un + (h / 2.0)*ku2, vn + (h / 2.0)*kv2);
		kv3 = func2(xn + h / 2.0, un + (h / 2.0)*ku2, vn + (h / 2.0)*kv2);

		ku4 = func1(xn + h, un + h * (ku1 - 2 * ku2 + 2 * ku3), vn + h * (kv1 - 2 * kv2 + 2 * kv3));
		kv4 = func2(xn + h, un + h * (ku1 - 2 * ku2 + 2 * ku3), vn + h * (kv1 - 2 * kv2 + 2 * kv3));

		un = un + h * (ku1 + 4.0 * ku3 + ku4) / 6.0;
		vn = vn + h * (kv1 + 4.0 * kv3 + kv4) / 6.0;

		return std::make_pair(un, vn);
	}
	std::pair<std::vector<double>, std::vector<double>> calculate()
	{
		arg.push_back(x0);
		auto result = std::make_pair(u0, v0);
		ures.push_back(u0);
		vres.push_back(v0);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(0.0);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		double xn = x0;
	//	double un = u0;
	//	double vn = v0;
		int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn <= xmax))
			{
				break;
			}
			else {
				if ((xn + h) > xmax)
				{
					while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
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
				ures.insert(ures.begin() + i + 1, result.first);
				vres.insert(vres.begin() + i + 1, result.second);
 				steps.insert(steps.begin() + i + 1, h);
				hinc.insert(hinc.begin() + i + 1, 0);
				hdec.insert(hdec.begin() + i + 1, 0);
				reswcap.insert(reswcap.begin() + i + 1, std::make_pair(0.0, 0.0));
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		N = i;
		Xn = arg[i];
		Un = ures[i];
		Vn = vres[i];
		if (i == n)
			hmin = 1;
		else hmin = i;
		hmax = 1;
		std::pair<std::vector<double>, std::vector<double>> res;
		res.first = ures;
		res.second = vres;
		return res;
	}
	std::pair<std::vector<double>, std::vector<double>> calculate_w_error()
	{
		auto result = std::make_pair(u0, v0);
		ss.push_back(0.0);
		arg.push_back(x0);
		ures.push_back(u0);
		vres.push_back(v0);
		steps.push_back(0.0);
//		steps.push_back(h);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		hinc.push_back(0);
		hdec.push_back(0);
		hinc.push_back(0);
		hdec.push_back(0);
		double xn = x0;
		double xhalf = x0;
		auto half = result;
		auto cap = result;
		auto next = result;
		double S, S1, S2;
		int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn <= xmax))
			{
				break;
			}
			else
			{
				next = RK4(xn, result.first, result.second, h);
				half = RK4(xn, result.first, result.second, h / 2.0);
				xhalf = xn + h / 2.0;
			    cap = RK4(xhalf, half.first, half.second, h / 2.0);
				S1 = abs(cap.first - next.first) / 15.0;
				S2 = abs(cap.second - next.second) / 15.0;
				S = std::max(S1, S2);
				if ((S >= (eps / 32.0)) && (S <= eps))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else
					{
						result = next;
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, dec);
					}

				}
				else if (S < (eps / 32.0))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						result = next;
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 2, ++inc);
						hdec.insert(hdec.begin() + i + 2, dec);
					}
				}
				else if (S > eps)
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						h = h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
				}
				i++;
				S *= 16.0;
				if (i == 1)
				{
					Smin = i;
					Smax = i;
					hmin = i;
					hmax = i;
				}
				
				else if (S < ss[Smin])
					Smin = i;
				else if (S > ss[Smax])
					Smax = i;
				if (h < steps[hmin])
					hmin = i;
				if (h > steps[hmax])
					hmax = i;
				reswcap.insert(reswcap.begin() + i, cap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				ures.insert(ures.begin() + i, result.first);
				vres.insert(vres.begin() + i, result.second);
			}
		}
		N = i;
		Xn = arg[i];
		Un = ures[i];
		Vn = vres[i];
		inc = hinc[i];
		dec = hdec[i];
		std::pair<std::vector<double>, std::vector<double>> res;
		res.first = ures;
		res.second = vres;
		return res;
	}
	friend std::ostream & operator<<(std::ostream &out, Pendulum &vc)
	{
		if (vc.ures.empty() || vc.ures.empty())
			out << "There are no calculated results yet.";
		else {
			out << "n  " << " h n-1  " << "    x    " << "    un    " << "    vn    " << "       S*      " << "inc " << "dec" << std::endl;
			for (int i = 0; i < vc.ures.size(); i++)
			{
				out << i << "  " << vc.steps[i] << "      " << vc.arg[i] << "    " << vc.ures[i] << "    " << vc.vres[i] << "    " << vc.ss[i] << "    " << vc.hinc[i] << "  " << vc.hdec[i] << std::endl;
			}
		}
		return out;
	}

};
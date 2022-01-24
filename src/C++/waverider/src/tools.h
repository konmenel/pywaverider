#pragma once
#include <iostream>
#include <vector>


using std::vector;


namespace wr 
{
	template <class T>
	bool not_active_condition(T x) { return false; };


	template <class in_type, class out_type>
	out_type runge_kutta(out_type(*func)(in_type, out_type), in_type x_bounds, out_type y0,
		double step = 1.0, bool (*condition)(in_type) = not_active_condition, uint32_t max_iter = 2000)
	{
		std::cout << func(x_bounds, y0) << step << condition(x_bounds) << max_iter << std::endl;
		return 0;

		/*
		out_type k1, k2, k3, k4;
		out_type y = y0;
		double x = x[0];

		// if first x larger than last x then we need to reverse condition `x < x_last` in while loop
		int corr = x[1] > x[0] ? 1 : -1;

		uint32_t n = 0;


		while (!condition<in_type>(y) && n < max_iter && corr * x < x[1])
		{
			k1 = step * func(x, y);
			k2 = step * func(x + 0.5 * step, y + 0.5 * k1);
			k3 = step * func(x + 0.5 * step, y + 0.5 * k2);
			k4 = step * func(x + step, y + k3);

			y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
			x += step;
			n++;
		}

		return y;
		*/
	}
}


typedef vector<double> out_t;
class RK45
{
private:
	double t0, tend, dt;
	out_t y0;
	vector<out_t(*)(double, out_t)> ode_func;
};
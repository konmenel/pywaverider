#pragma once
#include "tools.h"


double foo(int x, double y)
{
	return 1.0;
}


bool bar(int x)
{
	return 0;
}


int main()
{
	int x = 1;
	double y = 1.0;
	double res = wr::runge_kutta<int, double>(foo, x, y, 1.0, bar, 2000);
	bool m = wr::not_active_condition<int>(1);
	return 0;
}

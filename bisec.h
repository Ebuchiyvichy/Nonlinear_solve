#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"

local	localisation_func(int func, int nbrh)
{
	local	l;
	double	a, b;
	double	h;
	double	x;
	
	l.nbr = 0;
	boaders(func, &a, &b);
	h = fabs(b - a) / HBISEC;
	l.intervals = new double*[nbrh];	/*создаю лишние ячейки*/
	for (int i = 0; i <= nbrh; i++)
		l.intervals[i] = new double[2];
	for (int i = 0; i < nbrh; i++)
	{
		x = a + h * i;
		if (f(x, func) * f(x + h, func) <= 0)
		{
			l.intervals[l.nbr][0] = x;
			l.intervals[l.nbr][1] = x + h;
			l.nbr++;
		}
	}
	return l;
}

double	bisection(int func, local l, int n)
{
	double	x, a, b;

	a = l.intervals[n][0];
	b = l.intervals[n][1];
	x = (a + b) / 2;
	while (fabs(b - a) > EPS)
	{
		if (f(func, a) * f(func, x) <= 0)
			b = x;
		else
			a = x;
		x = (a + b) / 2;
	}
	return x;
}

double	D(int func, int x)
{
	return (f(x + EPS, func) - f(x, func)) / EPS;
}

double	Newton(int func, local l, int n)
{
	double	x1, x2, a, b;

	a = l.intervals[n][0];
	b = l.intervals[n][1];

	x2 = (f(a, func) * b - f(b, func) * a) / (f(a, func) - f(b, func));
	if (func == 3)
	{
		std::cout << "Trouble: ¯ \ _ (ツ) _ / ¯" << std::endl;
		return (x2);
	}
	do
	{
		x1 = x2;
		x2 = x1 - f(x1, func) / D(func, x1);
	} while (fabs(x1 - x2) > EPS);
	return (x2);
}

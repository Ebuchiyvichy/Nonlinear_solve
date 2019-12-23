#pragma once
#include <fstream>
#include <vector>
#include <cmath>

double			EPS = 10.e-6;
static int		HBISEC = 1;
int				TEST = 1;

typedef struct local
{
	int			nbr;
	double 		**intervals;
};

double	f(double x, int func)
{
	if (func == 1)
		return ((x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75));
	else if (func == 2)
		return sqrt(x + 1) - 1;
	else if (func == 3)
		return 35 * pow(x, 3) - 67 * pow(x, 2) - 3 * x + 3;
}

double	*f(double *x, int func)
{
	double	*y;

	if (func == -1)
	{
		y = new double[2];
		y[1] = x[1] * x[1] - x[2] * x[2] - 15;
		y[2] = x[1] * x[2] + 4;
		return y;
	}
	else if (func == -2)
	{
		y = new double[2];
		y[1] = x[1] * x[1] + x[2] * x[2] + x[1] + x[2] - 8;
		y[2] = x[1] * x[1] + x[2] * x[2] + x[1] * x[2] - 7;
		return y;
	}
}

void	boaders(int func, double *a, double *b)
{
	if (func == 1)
	{
		*a = 0.0;
		*b = 1.0;
	}
	else if (func == 2)
	{
		*a = -1.0;
		*b = 10.0;
	}
	else if (func == 3)
	{
		*a = 0.0;
		*b = 1.0;
	}
	else if (func == -1)
	{
		*a = 10;
		*b = 10;
	}
	else if (func == -2)
	{
		*a = 10;
		*b = 10;
	}
}
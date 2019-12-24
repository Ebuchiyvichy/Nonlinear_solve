#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"

double	D(int func, std::vector<double> x, int n, int m)
{
	std::vector<double>	tmp(NBR);

	for (int i = 0; i < NBR; i++)
		tmp[i] = x[i];
	tmp[m] = x[m] + EPS;
	return (f(tmp, func)[n] - f(x, func)[n]) / EPS;
}

double	D_cen(int func, std::vector<double> x, int n, int m)
{
	std::vector<double>	tmp1(NBR);
	std::vector<double>	tmp2(NBR);

	for (int i = 0; i < NBR; i++)
		tmp1[i] = x[i];
	tmp1[m] = x[m] + EPS;
	for (int i = 0; i < NBR; i++)
		tmp2[i] = x[i];
	tmp2[m] = x[m] - EPS;
	return (f(tmp1, func)[n] - f(tmp2, func)[n]) / (2 * EPS);
}

Matrix	Jacoby_matr(int func, std::vector<double> x)
{
	Matrix Jac(NBR);

	for (int i = 0; i < NBR; i++)
		for (int j = 0; j < NBR; j++)
			Jac.value[i][j] = D(func, x, i, j);
	return Jac;
}

std::vector<double>	Newton(int func, std::vector<double> x)
{
	std::vector<double>	xk(NBR);
	Matrix	J_(NBR);
	Matrix	R(NBR);
	Matrix	T(NBR);

	Matrix	J = Jacoby_matr(func, x);
	R = J;
	T.onebyone();
	T.QR_find_x(R);
	J_.inverse_matrix(R, T);
	do
	{
		xk = x;
		J = Jacoby_matr(func, xk);
		R = J;
		T.onebyone();
		T.QR_find_x(R);
		J_.inverse_matrix(R, T);
		x = diff(xk, multi_vect(f(xk, func), J_));
	} while (norm(xk, x) > EPS);
	return x;
}

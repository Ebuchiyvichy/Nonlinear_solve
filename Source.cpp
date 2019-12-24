#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "bisec.h"
#include "newton.h"

int	main()
{	
	if (TEST > 0)
	{
		local	l;
		int		tmp = 0;
		double	*x;

		l = localisation_func(TEST, HBISEC);
		while (l.nbr != tmp)		/*���������� ���������� ������ � ������ ���������� ��� ���*/
		{
			tmp = l.nbr;
			HBISEC += 3;
			l = localisation_func(TEST, HBISEC);
		}
		std::cout << "Number of solv is " << l.nbr << std::endl;

		x = new double[l.nbr];

		for (int i = 0; i < l.nbr; i++)
			x[i] = Bisection(TEST, l, i);	//��������� �� 0,04
		std::cout << "Bisection method:" << std::endl;
		for (int i = 0; i < l.nbr; i++)
			std::cout << "x" << i << " = " << x[i] << std::endl;

		for (int i = 0; i < l.nbr; i++)
			x[i] = Newton(TEST, l, i);		//��������� � 4� ����� � 1� �����
		std::cout << "Newton method:" << std::endl;
		for (int i = 0; i < l.nbr; i++)
			std::cout << "x" << i << " = " << x[i] << std::endl;

		for (int i = 0; i < l.nbr; i++)
			x[i] = Hords(TEST, l, i);		//��������� � 4� �����
		std::cout << "Hords method:" << std::endl;
		for (int i = 0; i < l.nbr; i++)
			std::cout << "x" << i << " = " << x[i] << std::endl;

		delete[] x;
	}
	else
	{
		Matrix	J(NBR);
		std::vector<double>	x(NBR);
		std::vector<double>	y(NBR);
		x[0] = 5;
		x[1] = 4;
	/*	J = Jacoby_matr(TEST, x);
		std::cout << "Matrix Jacoby:" << std::endl;
		J.print();*/
		y = Newton(TEST, x);
		for (int i = 0; i < NBR; i++)
			std::cout << y[i] << std::endl;
	}

	system("pause");
	return (0);
}
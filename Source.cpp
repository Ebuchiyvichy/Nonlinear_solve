#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "bisec.h"
#include "Reley.h"

int	main()
{	
	if (TEST > 0)
	{
		local	l;
		int		tmp = 0;
		double	*x;

		l = localisation_func(TEST, HBISEC);
/*		while (l.nbr != tmp)		//���������� ���������� ������ � ������ ���������� ��� ���
		{
			tmp = l.nbr;
			HBISEC += 3;
			l = localisation_func(TEST, HBISEC);
		}
*/		std::cout << "Number of solv is " << l.nbr << std::endl;

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
		std::vector<double>	a(N + 1);
		std::vector<double> b(N + 1);
		std::vector<double>	tmp1(NBR);
		std::vector<double>	tmp2(NBR);
		std::ofstream fout;
		int					h1, h2;

		x[0] = 5;	//��� ������� ����� - ������ (4 -1)
		x[1] = 4;	//��� ������� ����� ������ (1 2)
		y = Newton(TEST, x);
		for (int i = 0; i < NBR; i++)
			std::cout << y[i] << std::endl;
		if (TEST == -1)
		{
			x[0] = -5;	// ������ ���� ������ (-4 1)
			x[1] = -4;
			y = Newton(TEST, x);
			for (int i = 0; i < NBR; i++)
				std::cout << y[i] << std::endl;
		}
		else
		{
			x[0] = -5;	
			x[1] = -4;	//������ ���� ������ (-3 1)
			y = Newton(TEST, x);
			for (int i = 0; i < NBR; i++)
				std::cout << y[i] << std::endl;
			x[0] = -3;
			x[1] = -6;	/*������ ���� ������ (1 -3)*/
			y = Newton(TEST, x);
			for (int i = 0; i < NBR; i++)
				std::cout << y[i] << std::endl;
			x[0] = 3;
			x[1] = 6;	/*������ ���� ������ (2 1)*/
			y = Newton(TEST, x);
			for (int i = 0; i < NBR; i++)
				std::cout << y[i] << std::endl;
		}


		//������� �������
		a = uniform_mesh(TEST, N, 'a');		//������� �����
		b = uniform_mesh(TEST, N, 'b');
		std::cout << "Newton pool" << std::endl;
		fout.open("mass.txt");
		fout.trunc;
		fout << N << '\t' << 0 << '\t' << 0 << std::endl;
		for (int i = 0; i <= N; i++)
		{	
			iter = 0;
			for (int j = 0; j <= N; j++)
			{
				iter = 0;
				//������� �� ������ ������� �����
				x[0] = a[i];
				x[1] = b[j];
				y = Newton_reley(TEST, x);
				fout << iter << '\t';
				for (int i = 0; i < NBR; i++)
					fout << x[i] << '\t';					//�������� �������� � ���� x
				fout << std::endl;
			}
		}
		fout.close();
	}

	system("pause");
	return (0);
}
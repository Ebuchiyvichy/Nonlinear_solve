#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "bisec.h"

int	main()
{
	local	l;
	int		nbr  = 0;
	double	*x;
	
	l = localisation_func(TEST, HBISEC);
	while (l.nbr != nbr)		/*���������� ���������� ������ � ������ ���������� ��� ���*/
	{
		nbr = l.nbr;
		HBISEC += 3;
		l = localisation_func(TEST, HBISEC);	
	}
	std::cout << "Number of solv is " << l.nbr << std::endl;
	
	x = new double[l.nbr];

	for (int i = 0; i < l.nbr; i++)
		x[i] = bisection(TEST, l, i);	//��������� �� 0,04
	std::cout << "Bisection method:" << std::endl;
	for (int i = 0; i < l.nbr; i++)
		std::cout << "x" << i << " = " << x[i] << std::endl;

	for (int i = 0; i < l.nbr; i++)
		x[i] = Newton(TEST, l, i);		//��������� � 4� �����
	std::cout << "Newton method:" << std::endl;
	for (int i = 0; i < l.nbr; i++)
		std::cout << "x" << i << " = " << x[i] << std::endl;


	system("pause");
	return (0);
}
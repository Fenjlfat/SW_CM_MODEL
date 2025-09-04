#pragma once
#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>

double GNUPLOT(std::vector<double>& vectorVEL, std::vector<double>& vectorPRESS, std::vector<double>& vectorDENS, std::vector<double>& vectorENERG, std::vector<double>& Z);

void MODEL(std::vector<double>& input, std::vector<double>& output, double* KOEFFICIENTS, double* PARAM_MODELING, int FORMES);

struct PARAMETRS_MSS
{
	std::string direct_string = "XX";
	std::string vel_string = "V100";
	int time = 100000;
	int length = 100;  //количество €чеек или разбиений
	int direct = 1; //направление сжати€  1=XX; 2=YY; 3=ZZ
	int metall = 3; //1 = Al; 2 = Mg; 3 = CU

	double XX = 20.00e-3;  //длина моделируемой системы
	double VEL0 = 100.00e0; //скорость поршн€

	double VEPS = 1.00e5;
	double COEF_FORMS_PORE = 1.00e0; //SPHERE = 1.00

	double POROSITY0 = 0.45;
	double TEMPERATURE = 300.00e0;


	double alpha = 1000;   //коэффициент в€зкости
	double betta = 9.0;    //второй коэффициент в€зкости
	double timesss = 1.500e-6;  //до какого момента времени моделируем
	double correct_press = 0.05;
};

struct PARAMETRS
{
	int columnsFile = 10;
	int metall{ 2 }; //1=al; 2=cu

	double veps{ 1.0e5 };

	long double deltaTime = 7.00e-12;
	long double AMOL = 26.9815386e0;
	long double pi = 3.14159265358979323846264;
	long double kA = 20e0; // annihilation koefficient
	long double kB = 1.38e-23;
	long double Ha0 = 3 * 1e-3;
	long double Hb0 = 3 * 1e-3;
	long double DNA = 6.022045000e23; //Avogadro
};
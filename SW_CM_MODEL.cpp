// SW_CM_MODEL.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// SW_CM_TEST.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// SHOCK_WAVE_CM.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// SW_MSS_ANN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include "header.h"

using namespace::std;

//запись в файл
double FILES(double times, int length, vector<double>& vectorVEL, vector<double>& vectorPRESS, vector<double>& vectorDENS, vector<double>& vectorENERG, vector<double>& Z)
{
	ofstream EXIT;
	PARAMETRS_MSS PARAMETR;
	char tempChar[40];
	string tempSTRING;

	times = times * 1e12;

	sprintf_s(tempChar, "%.0f", times);

	tempSTRING = tempChar;

	EXIT.open("CUBE_CU_MSS_" + PARAMETR.direct_string + "_" + PARAMETR.vel_string + "_POR06_t" + tempSTRING + ".txt");
	EXIT << "NTS " << "LL(mkm) " << " STRESS(GPa) " << " " << " VEL(m / s) " << " DENS(kg/m3) " << " ENERGY(J) " << endl;
	for (size_t i = 0; i < length; i++)
	{
		EXIT << i << " " << Z[i]<< " " << vectorPRESS[i] * 1e-9 << " " << vectorVEL[i] << " " << vectorDENS[i] << " " << vectorENERG[i] << endl;

	}
	EXIT.close();
	return 0;
}

double minimum(double x, double y)
{
	if (x < y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

int main()
{
	PARAMETRS_MSS PAR_MSS;


	double TIME = 0e0, dx = 1e0, dt = 0e0, etta = 0e0;
	double GG = 0.00e0, KK = 0.00e0, Rho = 0e0;
	double DVEL = 0.00e0;
	double DEF = 0.00e0, DEF_XX = 0.00e0, DEF_YY = 0.00e0, DEF_ZZ = 0.00e0;
	string  times_file;

	//vector<vector<double>> VELOCITY(time, vector<double>(length + 1, 0));
	vector<double> VELOCITY0(PAR_MSS.length + 1, 0);
	vector<double> VELOCITY1(PAR_MSS.length + 1, 0);
	//vector<vector<double>> PRESS(time, vector<double>(length + 1, 0));
	vector<double> PRESS0(PAR_MSS.length + 1, 0);
	vector<double> PRESS1(PAR_MSS.length + 1, 0);
	//vector<vector<double>> DENSITY(time, vector<double>(length + 1, 0));
	vector<double> DENSITY0(PAR_MSS.length + 1, 0);
	vector<double> DENSITY1(PAR_MSS.length + 1, 0);
	//vector<vector<double>> STRESS(time, vector<double>(length + 1, 0));
	//vector<vector<double>> ENERGY(time, vector<double>(length + 1, 0));
	vector<double> ENERGY0(PAR_MSS.length + 1, 0);
	vector<double> ENERGY1(PAR_MSS.length + 1, 0);
	//vector<vector<double>> KSIG(time, vector<double>(length + 1, 0));
	vector<double> KSIG0(PAR_MSS.length + 1, 0);
	vector<double> KSIG1(PAR_MSS.length + 1, 0);
	vector<double> Z(PAR_MSS.length + 1, 0);//координаты
	//vector<vector<double>> h(time, vector<double>(length + 1, 0));//ширина ячеек
	vector<double> h0(PAR_MSS.length + 1, 0);
	vector<double> h1(PAR_MSS.length + 1, 0);
	vector<double> DEF11(PAR_MSS.length + 1, 0);
	vector<double> HH0(PAR_MSS.length + 1, 0);
	vector<double> mass(PAR_MSS.length + 1, 0);
	vector<double> ct(PAR_MSS.length + 1, 0);//скорость звука

	vector<double> input(32, 0.00);
	vector<double> output(32, 0.00);

	//AL
	if (PAR_MSS.metall == 1)
	{
		GG = 25.9e9 - 8.8e6 * PAR_MSS.TEMPERATURE; //AL
		KK = 79.9e9 + 1.9e7 * PAR_MSS.TEMPERATURE - 3.7e4 * pow(PAR_MSS.TEMPERATURE, 2); //AL
		Rho = 2700;
	}
	//MG
	if (PAR_MSS.metall == 2)
	{
		GG = 13.621745e9 - 1.032905e7 * PAR_MSS.TEMPERATURE; //MG
		KK = 35.186875e9 + 2.6505e7 * PAR_MSS.TEMPERATURE - 3.46375e-4 * pow(PAR_MSS.TEMPERATURE, 2); //MG
		Rho = 1738;
	}

	//CU переделать
	if (PAR_MSS.metall == 3)
	{
		GG = 40.7e9 - 8.8e6 * PAR_MSS.TEMPERATURE; //CU
		KK = 79.9e9 + 1.9e7 * PAR_MSS.TEMPERATURE - 3.7e4 * pow(PAR_MSS.TEMPERATURE, 2); //AL
		Rho = 8940;
	}

	dx = PAR_MSS.XX / PAR_MSS.length;
	Z[0] = 0.00e0;

	//распределение скорости по материалу
	for (int i = 0; i <= PAR_MSS.length; i++)
	{
		ct[i] = sqrt(GG / Rho);
		//Z[i] = Z[i-1]+1.00e-3/length;
		Z[i] = i * dx;
		if (i < 10)
		{
			//VELOCITY[0][i] = 300.00e0;
			VELOCITY0[i] = 0.00e0;
		}
		else
		{
			//VELOCITY[0][i] = 0.00e0;
			VELOCITY0[i] = 0.00e0;
		}
	}

	//распределение массы и плотности по ячейкам
	for (int i = 0; i < PAR_MSS.length; i++)
	{
		h0[i] = Z[i + 1] - Z[i];
		HH0[i] = h0[i];
		mass[i] = (h0[i] * Rho) * (1 - PAR_MSS.POROSITY0);
		DENSITY0[i] = mass[i] / h0[i];
	}
	std::cout << "in loop" << std::endl;
	int n = 0;
	while (TIME < PAR_MSS.timesss)
	{
		n = n + 1;
		//time step
		dt = 1.00e-2 * h0[0] / ct[0];
		for (int i = 0; i < PAR_MSS.length; i++)
		{
			dt = minimum(dt, 1.00e-2 * h0[i] / ct[i]);

		}
		//dt = 1e-10;
		TIME = TIME + dt;

		//viscosity
		for (int i = 0; i < PAR_MSS.length; i++)
		{

			//cout << PRESS[n-1][i] << "   ";
			DVEL = VELOCITY0[i + 1] - VELOCITY0[i];
			if (DVEL < 0.00e0)
			{
				//alpha = 300.00e0;  //ЧТО ЗА КОЭФФИЦИЕНТ
				etta = PAR_MSS.betta * h0[i] * Rho * PAR_MSS.alpha;
				KSIG0[i] = -etta * ((VELOCITY0[i + 1] - VELOCITY0[i]) / (h0[i]));
				//KSIG0[i] =0.00e0;
			}
			else
			{
				KSIG0[i] = 0.00e0;
			}

		}
		//вычисление скорости вещества с вязкостью
		for (int i = 0; i <= PAR_MSS.length; i++)
		{
			if (i == 0)
			{
				//VELOCITY1[i] = VELOCITY0[i] + dt * (-PRESS0[i] - KSIG0[i]) / (0.5 * (mass[i]));
				VELOCITY1[i] = PAR_MSS.VEL0;
			}
			else
			{
				if (i == PAR_MSS.length - 1)
					//VELOCITY1[i] = VELOCITY0[i] + dt * (-(-PRESS0[i] - KSIG0[i])) / (0.5 * (mass[i]));
					VELOCITY1[i] = VELOCITY0[i] + dt * (-PRESS0[i] - KSIG0[i] - (-PRESS0[i - 1] - KSIG0[i - 1])) / ((mass[i] + mass[i - 1]) * 0.5);
				else
					VELOCITY1[i] = VELOCITY0[i] + dt * (-PRESS0[i] - KSIG0[i] - (-PRESS0[i - 1] - KSIG0[i - 1])) / ((mass[i] + mass[i - 1]) * 0.5);
			}
			Z[i] = Z[i] + dt * (VELOCITY1[i] + VELOCITY0[i]) * 0.5;
			//cout << " i=" << i << "VEL = " << VELOCITY[n][i]<<endl;
		}

		//ширина ячейки
		for (int i = 0; i < PAR_MSS.length; i++)
		{
			h1[i] = Z[i + 1] - Z[i];
			DEF11[i] = -1 * (h1[i] - HH0[i]) / HH0[i];
		}

		//вычисление напряжения в ячейке
		for (int i = 0; i < PAR_MSS.length; i++)
		{
			//заполнение входного вектора
			input[0] = abs(DEF11[i]);
			input[1] = PAR_MSS.POROSITY0;
			input[2] = PAR_MSS.COEF_FORMS_PORE;
			input[3] = PAR_MSS.TEMPERATURE;
			input[4] = PAR_MSS.VEPS;

			double VEPS1 = PAR_MSS.VEPS;
			double VEPS2 = 0.00e0;
			double TEMPERATURE = PAR_MSS.TEMPERATURE;
			double KOEF_PORE;
			int FORMES = 2;
			if (FORMES == 1)
			{
				KOEF_PORE = 3.14159265358979323846264 / 6; //sphere
			}
			if (FORMES == 2)  //cube
			{
				KOEF_PORE = 1.00e0;
			}
			double XX = 1.2e-3;
			double ZZ = 1.2e-3;
			
			double EPSS = 0.501251;
			double EPSD = 1.06359;
			double GAMMA = 1.30874;
			double KANIH = 12.8492;
			double TAYLOR = 0.693304;
			
			
			double KOEFFICIENTS[5] = { EPSS, EPSD, GAMMA, KANIH, TAYLOR };
			double PARAM_MODELING[6] = { VEPS1, VEPS2, TEMPERATURE, KOEF_PORE, XX, ZZ };

			//вычисление модели и получение выходного вектора
			MODEL(input, output, KOEFFICIENTS, PARAM_MODELING, FORMES);

			PRESS1[i] = output[0] * -1e9 * (1 - output[2]) + output[0]*PAR_MSS.correct_press * 1e9;
			//PRESS1[i] = -KK * ((h1[i] - HH0[i]) / HH0[i]);
			//PRESS1[i] = PERCEPTRON(DEF, DEF_XX, DEF_YY, DEF_ZZ, COEF_FORMS_PORE, INIT_POROSITY, INITIAL_TEMP, VEPS) * 1e9;
		}
		//getchar();


		//вычисление плотности ячейки и скорсти звука в ней
		for (int i = 0; i < PAR_MSS.length; i++)
		{
			DENSITY1[i] = mass[i] / h1[i];
			ct[i] = sqrt(KK / DENSITY1[i]);
			//cout << DENSITY[n][i] << "   ";
		}

		//энергия
		for (int i = 0; i <= PAR_MSS.length - 1; i++)
		{
			ENERGY1[i] = ENERGY0[i] - (PRESS0[i] / pow(DENSITY0[i], 2)) * (DENSITY1[i] - DENSITY0[i]);
		}
		//перезапись текущих значений
		for (int i = 0; i < PAR_MSS.length; i++)
		{
			PRESS0[i] = PRESS1[i];
			VELOCITY0[i] = VELOCITY1[i];
			DENSITY0[i] = DENSITY1[i];
			ENERGY0[i] = ENERGY1[i];
			h0[i] = h1[i];
			KSIG0[i] = KSIG1[i];
		}

		cout << endl << " t=" << TIME << "	n=" << n << endl;

		//вывод графиков для определенного времени
		if (TIME > 8.00e-7 && TIME < 8.011e-7)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 8.50e-7 && TIME < 8.511e-7)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 9.00e-7 && TIME < 9.011e-7)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 9.50e-7 && TIME < 9.511e-7)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 1.00e-6 && TIME < 1.011e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 1.10e-6 && TIME < 1.111e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 1.5e-6 && TIME < 1.514e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 2.00e-6 && TIME < 2.011e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 2.10e-6 && TIME < 2.11e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 2.50e-6 && TIME < 2.511e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 3.00e-6 && TIME < 3.01e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 3.10e-6 && TIME < 3.11e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 3.40e-6 && TIME < 3.41e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		if (TIME > 3.50e-6 && TIME < 3.511e-6)
		{
			double tttt = TIME;
			FILES(tttt, PAR_MSS.length, VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
			//GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
		}
		//getchar();
	}
	//getchar();
	GNUPLOT(VELOCITY1, PRESS1, DENSITY1, ENERGY1, Z);
	std::cout << "Hello World!\n";
}

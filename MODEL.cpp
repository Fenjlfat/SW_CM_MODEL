#include "header.h"

double MINIMUM(double LA, double LHA)
{
	if (LA < LHA)
	{
		return LA;
	}
	else
	{
		return LHA;
	}
}
double MAXIMUM(double LA, double LHA)
{
	if (LA > LHA)
	{
		return LA;
	}
	else
	{
		return LHA;
	}
}
double SIGNUM(double x)
{

	if (x < 0)
	{
		return -1;
	}
	else
	{
		return 1;
	}

}

void MODEL(std::vector<double>& input, std::vector<double>& output, double* KOEFFICIENTS, double* PARAM_MODELING, int FORMES)
{
	PARAMETRS PAR;
	double TIME = 0.00e0, t = 0.00e0, dt = PAR.deltaTime;
	double eps1 = 0.e0, eps2 = 0.e0;
	double LL0 = 0.00e0;
	double LLa = 0.e0, LLb = 0.e0;
	double L11 = 0.e0, L22 = 0.e0;
	double Fp11 = 0.e0, Fp22 = 0.e0;
	double Fd11 = 0.e0, Fd22 = 0.e0;
	double AA1 = 0.e0, AA2 = 0.e0;
	double LEN1 = 0.e0, LEN2 = 0.e0;
	double AC1 = 0.e0, AC2 = 0.e0;
	double VOLP0 = 0.e0, VOLP = 0.e0;
	double VOLS0 = 0.e0, VOLS = 0.0e0;
	double VOLC = 0.e0, ACR = 0.0e0;
	double RHOD = 0.e0, ALF = 0.e0;
	double vFd11_0 = 0.e0, vFd22_0 = 0.e0;
	double F11 = 0.e0, F22 = 0.e0;
	double Fe11 = 0.e0, Fe22 = 0.e0;
	double PRES = 0.e0, SDEV = 0.e0;
	double SIG1 = 0.e0, SIG2 = 0.e0;
	double PF1 = 0.e0, PF2 = 0.e0;
	double YY = 0.e0, ETA1 = 0.e0;
	double TAU1 = 0.e0, TAU2 = 0.e0;
	double Wpl1 = 0.e0, Wpl2 = 0.e0, Wpl = 0.e0;
	double AA1_0 = 0.e0, AA2_0 = 0.e0;
	double AC1_0 = 0.e0, AC2_0 = 0.e0;
	double SIG_NUCL = 0.0e0, PROB_NUCL = 0.e0;

	double C_FCC = 0.00e0, Ct = 0.00e0, AMOL = 0.00e0;
	double m1 = 0.00e0, d1 = 0.00e0;
	double DNC = 0.00e0, BURG = 0.00e0;
	double GG = 0.00e0, KK = 0.00e0;
	double DENS = 0.00e0, BTR = 0.00e0;


	double VEPS1 = PARAM_MODELING[0];
	double VEPS2 = PARAM_MODELING[1];
	double TEMPERATURE = PARAM_MODELING[2];
	double KOEF_PORE = PARAM_MODELING[3];
	double XX = PARAM_MODELING[4];
	double ZZ = PARAM_MODELING[5];

	double EPSS = KOEFFICIENTS[0];
	double EPSD = KOEFFICIENTS[1];
	double GAMMA = KOEFFICIENTS[2];
	double KANIH = KOEFFICIENTS[3];
	double TAYLOR = KOEFFICIENTS[4];

	//GAM = gamma1;  //!*(1.e0 + 4.13845369e-4 * TEMP);

	//elastic constants and density
	if (PAR.metall == 1) //AL
	{
		AMOL = 26.9815386e0;
		GG = 25.9e9 - 8.8e6 * TEMPERATURE;
		KK = 79.9e9 + 1.9e7 * TEMPERATURE - 3.7e4 * pow(TEMPERATURE, 2);
		DENS = 2.7e3 - 8.7e-2 * TEMPERATURE - 5.41e-5 * pow(TEMPERATURE, 2);
		BTR = 0.674e-5 + 2.55e-8 * TEMPERATURE;
		EPSD = EPSD * (1.e0 - 3.4e-4 * TEMPERATURE);
		EPSS = EPSS * (1.e0 - 3.4e-4 * TEMPERATURE);
	}
	if (PAR.metall == 2) //CU
	{
		AMOL = 63.546e0;
		GG = 49.79957143e9 - 2.150571429e7 * TEMPERATURE;
		KK = 145.5314643e9 - 1.24547619e7 * TEMPERATURE - 1.963690476e4 * pow(TEMPERATURE, 2);
		DENS = 8.890885714e3 - 3.575e-1 * TEMPERATURE - 1.192857143e-4 * pow(TEMPERATURE, 2);
		BTR = 1.2e-5;
		EPSD = EPSD * (1.e0 - 4.13845369e-4 * TEMPERATURE);
		EPSS = EPSS * (1.e0 - 3.4e-4 * TEMPERATURE);
	}

	C_FCC = 0.5e0 * sqrt(2.e0) * pow((4.e0), (0.333333));
	Ct = sqrt(GG / DENS);

	//Burgers vector
	m1 = AMOL * 1.e-3 / PAR.DNA;
	DNC = (DENS / m1);
	d1 = pow(DNC, 0.333333);
	BURG = C_FCC / d1;



	TIME = 0.e0;
	eps1 = 0.e0;
	eps2 = 0.e0;
	LL0 = PAR.Ha0;

	Fp11 = 1.e0;
	Fp22 = 1.e0;

	Fd11 = 1.e0;
	Fd22 = 1.e0;

	AA1 = XX;
	AA2 = ZZ;

	LEN1 = MINIMUM(AA1, LLa - AA1);
	LEN2 = MINIMUM(AA2, LLb - AA2);
	//!TAU1 = (-GAM) / LEN1;
	//!TAU2 = (-GAM) / LEN2;
	//!AC1 = AA1 * (1.d0 + TAU1 / (2.d0 * GG));
	//!AC2 = AA2 * (1.d0 + TAU2 / (2.d0 * GG));

	AC1 = AA1 - GAMMA / (2.e0 * GG);
	AC2 = AA2 - GAMMA / (2.e0 * GG);


	VOLP = KOEF_PORE * AA1 * AA2 * AA2;
	VOLC = KOEF_PORE * AC1 * AC2 * AC2;

	PROB_NUCL = 0.e0;
	RHOD = 0.e16;
	ALF = VOLC / (pow(LL0, 3));
	VOLS0 = pow(LL0, 3) - VOLP;

	vFd11_0 = 0.e0;
	vFd22_0 = 0.e0;

	VOLS = VOLS0;

	t = 0.e0;

	int j = 0;

	double epsM = input[0];

	while (eps1 + 2 * eps2 < epsM)
	{
		TIME = TIME + dt;

		eps1 = eps1 + VEPS1 * dt;
		eps2 = eps2 + VEPS2 * dt;
		LLa = LL0 * (1.e0 - eps1);
		LLb = LL0 * (1.e0 - eps2);

		F11 = (1.e0 - eps1);
		F22 = (1.e0 - eps2);
		Fe11 = F11 / (Fp11 * Fd11);
		Fe22 = F22 / (Fp22 * Fd11);
		L11 = 0.5e0 * (pow(Fe11, 2) - 1.e0);
		L22 = 0.5e0 * (pow(Fe22, 2) - 1.e0);

		PRES = -KK * (L11 + 2.e0 * L22);
		SDEV = (4.e0 / 3.e0) * GG * (L11 - L22);
		SIG1 = -PRES + SDEV;
		SIG2 = -PRES - SDEV * 0.5e0;

		LEN1 = MINIMUM(AA1, LLa - AA1);
		LEN2 = MINIMUM(AA2, LLb - AA2);

		//Burgers vector
		m1 = AMOL * 1.e-3 / PAR.DNA;
		DNC = (DENS / m1);
		d1 = pow(DNC, 0.333333);
		BURG = C_FCC / d1;
		Ct = sqrt(GG / DENS);

		if (FORMES == 1) //sphera 
		{
			PF1 = PAR.pi * AA2;
			PF2 = PAR.pi * AA1;
		}

		if (FORMES == 2) //cube
		{
			PF1 = 4.e0 * AA2;
			PF2 = 2.e0 * (AA1 + AA2);
		}
		//collapse of pores
		YY = 30.e6 + TAYLOR * GG * BURG * sqrt(RHOD);//CU

		ETA1 = pow(BURG, 2) * RHOD / (4.e0 * BTR);

		TAU1 = 3.e0 * GG * (1.e0 + PRES * (1.e0 - ALF) / (3.e0 * KK) - 2.e0 * GG * AA1 / (2.e0 * GG * AA1 + GAMMA) * (1.e0 + PRES * (1.e0 / (4.e0 * GG) + 1.e0 / (3.e0 * KK))));
		TAU2 = 3.e0 * GG * (1.e0 + PRES * (1.e0 - ALF) / (3.e0 * KK) - 2.e0 * GG * AA2 / (2.e0 * GG * AA2 + GAMMA) * (1.e0 + PRES * (1.e0 / (4.e0 * GG) + 1.e0 / (3.e0 * KK))));

		if (abs(TAU1) > 0.5e0 * YY)
		{
			Wpl1 = -ETA1 * abs(TAU1 + (0.75e0 * SDEV) * SIGNUM(TAU1) - 0.5e0 * YY * SIGNUM(TAU1));
		}
		else
		{
			Wpl1 = 0.e0;
			//Wpl1 = ETA1 * abs(TAU1 + (0.75e0 * SDEV) * SIGNUM(TAU1) - 0.5e0 * YY * SIGNUM(TAU1));
		}

		if (abs(TAU2) > 0.5e0 * YY)
		{
			Wpl2 = -ETA1 * abs(TAU2 + (0.75e0 * SDEV) * SIGNUM(TAU2) - 0.5e0 * YY * SIGNUM(TAU2));
		}
		else
		{
			Wpl2 = 0.e0;
			//Wpl2 = ETA1 * abs(TAU2 + (0.75e0 * SDEV) * SIGNUM(TAU2) - 0.5e0 * YY * SIGNUM(TAU2));
		}

		AA1_0 = AA1;
		AA2_0 = AA2;

		AA1 = AA1 * (1.e0 + dt * Wpl1);
		AA2 = AA2 * (1.e0 + dt * Wpl2);

		AC1_0 = AC1;
		AC2_0 = AC2;

		AC1 = (AA1 - 0.5e0 * GAMMA / GG) / (1.e0 + PRES * (1.e0 - ALF) * (1.e0 / (4.e0 * GG) + 1.e0 / (3.e0 * KK)));
		AC2 = (AA2 - 0.5e0 * GAMMA / GG) / (1.e0 + PRES * (1.e0 - ALF) * (1.e0 / (4.e0 * GG) + 1.e0 / (3.e0 * KK)));

		VOLP0 = VOLP;
		VOLP = KOEF_PORE * AA1 * AA2 * AA2;
		VOLC = KOEF_PORE * AC1 * AC2 * AC2;
		ALF = VOLC / (LLa * LLb * LLb);
		VOLS = LLa * LLb * LLb - VOLP;
		//VOLS = LLa * LLb * LLb;
		/*!Fd11 = Fd11 * (1.d0 + KOEF_FORM * (AA1 - AA1_0) * AA2 * *2 / (VOLS0 * Fe11 * Fp11 * (Fe22 * Fp22) * *2)) + &
		!(1.d0 - Fd11 * Fd22 * *2 * LL0 * *3 / VOLS0) * (-veps1 * dt) / (Fe11 * Fp11)
		!Fd22 = Fd22 * (1.d0 + KOEF_FORM * (AA2 - AA2_0) * AA1 * AA2 / (VOLS0 * Fe11 * Fp11 * (Fe22 * Fp22) * *2)) + &
		!(1.d0 - Fd11 * Fd22 * *2 * LL0 * *3 / VOLS0) * (-veps2 * dt) / (Fe22 * Fp22)
		*/

		Fd11 = Fd11 * (1.e0 + KOEF_PORE * (0.5e0 * (AC1 - AC1_0) + 0.5e0 * (AA1 - AA1_0)) * AA2 * AA2
			/ (VOLS0 * Fe11 * Fp11 * pow((Fe22 * Fp22), 2))) + (1.e0 - Fd11 * pow(Fd22, 2) * pow(LL0, 3) / VOLS0) * (-VEPS1 * dt) / (Fe11 * Fp11);
		Fd22 = Fd22 * (1.e0 + KOEF_PORE * (0.5e0 * (AC2 - AC2_0) + 0.5e0 * (AA2 - AA2_0)) * AA1 * AA2
			/ (VOLS0 * Fe11 * Fp11 * pow((Fe22 * Fp22), 2))) + (1.e0 - Fd11 * pow(Fd22, 2) * pow(LL0, 3) / VOLS0) * (-VEPS2 * dt) / (Fe22 * Fp22);

		/*
		!Fd11 = Fd11 * (1.d0 + KOEF_FORM * ((AC1 - AC1_0)) * AA2 * *2 / (VOLS0 * Fe11 * Fp11 * (Fe22 * Fp22) * *2)) + &
		!(1.d0 - Fd11 * Fd22 * *2 * LL0 * *3 / VOLS0) * (-veps1 * dt) / (Fe11 * Fp11)
		!Fd22 = Fd22 * (1.d0 + KOEF_FORM * ((AC2 - AC2_0)) * AA1 * AA2 / (VOLS0 * Fe11 * Fp11 * (Fe22 * Fp22) * *2)) + &
		!(1.d0 - Fd11 * Fd22 * *2 * LL0 * *3 / VOLS0) * (-veps2 * dt) / (Fe22 * Fp22)
		*/

		//plastic deformation
		if (abs(SDEV) > 0.5e0 * YY * 4.e0 / 3.e0)
		{
			Wpl = ETA1 * (SDEV - 0.5e0 * YY * SIGNUM(SDEV) * 4.e0 / 3.e0);
		}
		else
		{
			Wpl = 0.e0;
		}

		Fp11 = Fp11 * exp(Wpl * dt);
		Fp22 = Fp22 * exp(-0.5e0 * Wpl * dt);

		//emission of dislocations
		SIG_NUCL = abs(0.75e0 * SDEV) + MAXIMUM(abs(TAU1), abs(TAU2));
		ACR = (EPSS * 1.6e-19 / BURG) / (BURG * SIG_NUCL);
		PROB_NUCL = PROB_NUCL + (2.e0 * PAR.pi * AA1 * AA2 * BURG * Ct / pow(ACR, 4)) * exp(-PAR.pi * (EPSS * 1.6e-19 / BURG) * ACR / (2.e0 * PAR.kB * TEMPERATURE));
		if (PROB_NUCL > 1.e0)
		{
			RHOD = RHOD + PAR.pi * ACR / (LLa * LLb * LLb); // !*PROB_NUCL
			PROB_NUCL = 0.e0;
		}

		RHOD = RHOD + dt * ((abs(Wpl * SDEV * 1.5e0) + (abs(Wpl1 * TAU1) + abs(Wpl2 * TAU2)) * ALF * 0.5e0) / (EPSD * 1.6e-19 / BURG) - KANIH * RHOD * (abs(Wpl) + (abs(Wpl1) + abs(Wpl2)) * ALF * 0.5e0));

		t = t + dt;
		output[0] = SIG1 / 1e9;
		output[1] = SIG2 / 1e9;
		output[2] = VOLP / VOLS;
	}
	j = 0;
	//getchar();
}
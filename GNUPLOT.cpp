#include "header.h"
//рисование графиков
double GNUPLOT(std::vector<double>& vectorVEL, std::vector<double>& vectorPRESS, std::vector<double>& vectorDENS, std::vector<double>& vectorENERG, std::vector<double>& Z)
{
	FILE* pipe = _popen("C:\\gnuplot\\bin\\gnuplot.exe", "w");
	if (pipe != NULL)
	{
		//PARAMETRS_MSS PARAMETR;
		std::ofstream DATAGNUPLOT;
		char tempChar[40];

		std::string tempSTRING;

		tempSTRING = tempChar;

		DATAGNUPLOT.open("DataGnuplot.txt");
		for (int i = 0; i < vectorPRESS.size() - 1; i++)
		{
			DATAGNUPLOT << i << " " << Z[i] << " " << vectorPRESS[i] * 1e-9 << " " << -vectorVEL[i] << " " << vectorDENS[i] << " " << vectorENERG[i] << std::endl;   //SIG11
		}
		DATAGNUPLOT.close();

		//fprintf(pipe, "plot  'E:\\WORK\\SPH\\PROGRAMMS\\C++\\SW_CM_MODEL\\SW_CM_MODEL\\DataGnuplot.txt' using 1:3 title 'MODEL' with lines linecolor 5\n");
		//fprintf(pipe, "plot 'E:\\WORK\\SPH\\SW_CM_TEST\\SW_CM_TEST\\DataGnuplot.txt' using 1:3 title 'MODEL' with lines linecolor 1,'E:\\WORK\\SPH\\SW_CM_TEST\\SW_CM_TEST\\400ms_2mks_2.txt' using 8:4 title 'SPH' with lines linecolor 5\n");
		fprintf(pipe, "plot 'E:\\WORK\\SPH\\PROGRAMMS\\C++\\SW_CM_MODEL\\SW_CM_MODEL\\DataGnuplot.txt' using 1:3 title 'MODEL' with lines linecolor 1,'E:\\WORK\\SPH\\PROGRAMMS\\C++\\SW_CM_MODEL\\SW_CM_MODEL\\400MS_25MKS.txt' using 8:4 title 'SPH' with lines linecolor 5\n");

		fflush(pipe);

		// ожидание нажатия клавиши
		std::cin.clear();
		std::cin.ignore(std::cin.rdbuf()->in_avail());
		std::cin.get();

#ifdef WIN32
		_pclose(pipe);
#else
		_pclose(pipe);
#endif
	}
	else puts("Could not open the file\n");
	//getchar();
	_pclose(pipe);
	return 0;
}
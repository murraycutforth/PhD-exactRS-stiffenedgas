#include "exact_RS_stiffenedgas.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <blitz/array.h>



double stiffenedgas_e (double rho, double p, double gamma, double pinf)
{
	return (p+gamma*pinf)/(rho*(gamma-1));
}



int main()
{
	// Output solution for each of the five Toro test cases

	for (int TC = 1; TC <= 5; TC++)
	{
		
		blitz::Array<double,1> WL (3);
		blitz::Array<double,1> WR (3);
		double t;
		std::string filename;
	
		double gammaL = 1.4;
		double gammaR = 1.4;
		double pinf_L = 0.0;
		double pinf_R = 0.0;
	
		if (TC == 1)
		{
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 0.1;
			t = 0.25;
			filename = "TTC1.dat";
		}
		else if (TC ==2)
		{
			WL(0) = 1.0;
			WL(1) = -2.0;
			WL(2) = 0.4;
			WR(0) = 1.0;
			WR(1) = 2.0;
			WR(2) = 0.4;
			t = 0.15;
			filename = "TTC2.dat";
		}
		else if (TC ==3)
		{
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1000.0;
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 0.01;
			t = 0.012;
			filename = "TTC3.dat";
		}
		else if (TC ==4)
		{
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 0.01;
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 100.0;
			t = 0.035;
			filename = "TTC4.dat";
		}
		else if (TC ==5)
		{
			WL(0) = 5.99924;
			WL(1) = 19.5975;
			WL(2) = 460.894;
			WR(0) = 5.99242;
			WR(1) = -6.19633;
			WR(2) = 46.0950;
			t = 0.035;
			filename = "TTC5.dat";
		}
	
	
		exact_rs_stiffenedgas RS (gammaL, gammaR, pinf_L, pinf_R);
		RS.solve_RP(WL,WR);
	
		std::cout << "Star state pressure calculated as " << RS.P_STAR << std::endl;
		std::cout << "Star state velocity calculated as " << RS.S_STAR << std::endl;
		std::cout << "Left star state density calculated as " << RS.rho_star_L << std::endl;
		std::cout << "Right state state density calculated as " << RS.rho_star_R << std::endl;
		std::cout << "Left shock speed calculated as " << RS.S_L << std::endl;
		std::cout << "Right shock speed calculated as " << RS.S_R << std::endl;
		std::cout << "Left rarefaction head speed calculated as " << RS.S_HL << std::endl;
		std::cout << "Left rarefaction tail speed calculated as " << RS.S_TL << std::endl;
		std::cout << "Right rarefaction head speed calculated as " << RS.S_HR << std::endl;
		std::cout << "Right rarefaction tail speed calculated as " << RS.S_TR << std::endl;
	
		double xmin = 0.0;
		double xmax = 1.0;
		double offset = 0.5;
		int numsamples = 10000;
		double delx = (xmax - xmin)/numsamples;
	
		std::ofstream outfile;
		outfile.open(filename);
	
		double x = xmin;
		while (x <= xmax)
		{
			double S = x/t;
			blitz::Array<double,1> soln (RS.sample_solution(WL, WR, S - offset/t));
			double thisgamma = S - offset/t < RS.S_STAR ? gammaL : gammaR;
			double thispinf = S - offset/t < RS.S_STAR ? pinf_L : pinf_R;
			outfile << x << " " << soln(0) << " " << soln(1) << " " << soln(2) << " " << stiffenedgas_e(soln(0), soln(2), thisgamma, thispinf) << std::endl;
			x += delx;
		}
	}

	return 0.0;
}

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

	for (int TC = 16; TC <= 17; TC++)
	{
		
		blitz::Array<double,1> WL (3);
		blitz::Array<double,1> WR (3);
		double t;
		std::string filename;
	
		double gammaL = 1.4;
		double gammaR = 1.4;
		double pinf_L = 0.0;
		double pinf_R = 0.0;
		double offset = 0.5;
		double xmin = 0.0;
		double xmax = 1.0;
	
		if (TC == 1)
		{
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 0.1;
			t = 0.25;
			filename = "./output/TTC1.dat";
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
			filename = "./output/TTC2.dat";
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
			filename = "./output/TTC3.dat";
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
			filename = "./output/TTC4.dat";
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
			filename = "./output/TTC5.dat";
		}
		else if (TC == 6)
		{
			// Air - Helium shock tube from Sambasivan 2009
			
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 0.1;
			t = 0.25;
			filename = "./output/Samb1.dat";
			gammaR = 1.667;
		}
		else if (TC == 7)
		{
			// Gaseous shock tube from rGFM
			
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 100000.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 10000.0;
			t = 0.0007;
			filename = "./output/NE1.dat";
			gammaR = 1.2;
		}
		else if (TC == 8)
		{
			// Air - water shock from rGFM
			
			WL(0) = 0.00596521;
			WL(1) = 911.8821;
			WL(2) = 1000.0;
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 1.0;
			gammaR = 7.15;
			pinf_R = 3309.0;
			t = 0.0007;
			filename = "./output/rGFM2.dat";
		}
		else if (TC == 9)
		{
			// Air - water jet from rGFM
			
			WL(0) = 1.0;
			WL(1) = 90.0;
			WL(2) = 1.0;
			WR(0) = 1000.0;
			WR(1) = 0.0;
			WR(2) = 1.0;
			gammaR = 7.15;
			pinf_R = 3309.0;
			t = 0.015;
			filename = "./output/rGFM4.dat";
			offset = 0.6;
		}
		else if (TC == 10)
		{
			// Reversed version of TC 9
			
			WR(0) = 1.0;
			WR(1) = -90.0;
			WR(2) = 1.0;
			WL(0) = 1000.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			gammaL = 7.15;
			pinf_L = 3309.0;
			t = 0.015;
			filename = "./output/rGFM_reversed.dat";
			offset = 0.4;
		}
		else if (TC == 11)
		{
			// Water - air shock from Saurel 1999
			
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 4.0;
			WR(0) = 0.05;
			WR(1) = 0.0;
			WR(2) = 0.0004;
			gammaL = 4.4;
			pinf_L = 2.4;
			t = 0.12;
			filename = "./output/Saurel1.dat";
			offset = 0.7;
		}
		else if (TC == 12)
		{
			// Reversed version of TC 11
			
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 4.0;
			WL(0) = 0.05;
			WL(1) = 0.0;
			WL(2) = 0.0004;
			gammaR = 4.4;
			pinf_R = 2.4;
			t = 0.12;
			filename = "./output/Saurel1_reversed.dat";
			offset = 0.3;
		}
		else if (TC == 13)
		{
			// Numerical experiment 2
			
			WL(0) = 3.175962;
			WL(1) = 9.434992;
			WL(2) = 100.0;
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 1.0;
			gammaL = 1.667;
			gammaR = 1.2;
			t = 0.045;
			filename = "./output/NE2.dat";
			offset = 0.2;
		}
		else if (TC == 14)
		{
			// Numerical experiment 3
			
			WL(0) = 0.00596521;
			WL(1) = 911.8821;
			WL(2) = 100.0;
			WR(0) = 1.0;
			WR(1) = 0.0;
			WR(2) = 1.0;
			gammaL = 1.4;
			gammaR = 7.15;
			pinf_L = 0.0;
			pinf_R = 3309.0;
			t = 0.0007;
			filename = "./output/NE3.dat";
			offset = 0.5;
		}
		else if (TC == 15)
		{
			// Gaseous shock tube from So/Hu/Adams 2012
			
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 0.1;
			t = 0.15;
			filename = "./output/ST3.dat";
			gammaR = 1.667;
		}
		else if (TC == 16)
		{
			// Gaseous shock tube from Garrick/Owkes/Regele 2016
			
			WL(0) = 1.0;
			WL(1) = 0.0;
			WL(2) = 1.0;
			WR(0) = 0.125;
			WR(1) = 0.0;
			WR(2) = 0.1;
			t = 0.14;
			filename = "./output/ST1.dat";
			gammaR = 2.4;
		}
		else if (TC == 17)
		{
			// Water - air shock tube from Murrone/Guillard 2004
			
			WL(0) = 1000.0;
			WL(1) = 0.0;
			WL(2) = 1000000000.0;
			WR(0) = 50.0;
			WR(1) = 0.0;
			WR(2) = 100000.0;
			t = 0.0009;
			filename = "./output/ST2.dat";
			offset = 0.7;
			gammaL = 4.4;
			gammaR = 1.4;
			pinf_L = 600000000.0;
			pinf_R = 0.0;
			xmin = -2.0;
			xmax = 2.0;
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
		std::cout << "Right rarefaction tail speed calculated as " << RS.S_TR << std::endl << std::endl;
	
		
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
			double thisz = S - offset/t < RS.S_STAR ? 1.0 : 0.0;
			outfile << x << " " << soln(0) << " " << soln(1) << " " << soln(2) << " " << stiffenedgas_e(soln(0), soln(2), thisgamma, thispinf) << " " << thisz << std::endl;
			x += delx;
		}
	}

	return 0.0;
}

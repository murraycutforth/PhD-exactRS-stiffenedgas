#include "exact_RS_stiffenedgas.hpp"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>



exact_rs_stiffenedgas :: exact_rs_stiffenedgas (double gamma_L, double gamma_R, double pinf_L, double pinf_R)
:
	gamma_L	(gamma_L),
	gamma_R	(gamma_R),
	pinf_L (pinf_L),
	pinf_R (pinf_R),
	S_STAR	(0.0),
	P_STAR	(0.0),
	rho_star_L	(0.0),
	rho_star_R	(0.0),
	S_L	(0.0),
	S_R	(0.0),
	S_HL	(0.0),
	S_TL	(0.0),
	S_HR	(0.0),
	S_TR	(0.0)
{}




void exact_rs_stiffenedgas :: solve_RP (const Eigen::Vector3d& W_L, const Eigen::Vector3d& W_R)
{
	assert(W_L(0) >= 0.0);
	assert(W_R(0) >= 0.0);
	// assert(W_L(2) >= 0.0); Since stiffened gases will often exhibit p<0..
	// assert(W_R(2) >= 0.0);

	
	// Calculate p_star

	P_STAR = find_p_star_newtonraphson(W_L(0), W_L(1), W_L(2), W_R(0), W_R(1), W_R(2));	
		
	
	// Calculate u_star

	S_STAR = 0.5*(W_L(1)+W_R(1)) + 0.5*(f(P_STAR,W_R(0),W_R(2),gamma_R,pinf_R) - f(P_STAR,W_L(0),W_L(2),gamma_L,pinf_L));


	// Solution now depends on character of 1st and 3rd waves

	if (P_STAR > W_L(2))
	{
		// Left shock
		
		rho_star_L = W_L(0)*((2.0*gamma_L*pinf_L + (gamma_L+1.0)*P_STAR + (gamma_L-1.0)*W_L(2))/(2.0*(W_L(2) + gamma_L*pinf_L) + (gamma_L-1.0)*P_STAR + (gamma_L-1.0)*W_L(2)));
		S_L = W_L(1) - (Q_K(P_STAR,W_L(0),W_L(2),gamma_L,pinf_L)/W_L(0));

	}
	else
	{
		// Left rarefaction

		rho_star_L = W_L(0)*std::pow((P_STAR + pinf_L)/(W_L(2) + pinf_L), 1.0/gamma_L);

		double a_L = a(W_L(0), W_L(2), gamma_L, pinf_L);
		double a_star_L = a_L*std::pow((P_STAR + pinf_L)/(W_L(2) + pinf_L), (gamma_L-1.0)/(2.0*gamma_L));

		S_HL = W_L(1) - a_L;
		S_TL = S_STAR - a_star_L;
	}

	if (P_STAR > W_R(2))
	{
		// Right shock
				
		rho_star_R = W_R(0)*((2.0*gamma_R*pinf_R + (gamma_R+1.0)*P_STAR + (gamma_R-1.0)*W_R(2))/(2.0*(W_R(2) + gamma_R*pinf_R) + (gamma_R-1.0)*P_STAR + (gamma_R-1.0)*W_R(2)));

		S_R = W_R(1) + (Q_K(P_STAR,W_R(0),W_R(2),gamma_R,pinf_R)/W_R(0));
	}
	else
	{
		// Right rarefaction

		rho_star_R = W_R(0)*std::pow((P_STAR + pinf_R)/(W_R(2) + pinf_R), 1.0/gamma_R);

		double a_R = a(W_R(0),W_R(2),gamma_R, pinf_R);
		double a_star_R = a_R*std::pow((P_STAR + pinf_R)/(W_R(2) + pinf_R), (gamma_R-1.0)/(2.0*gamma_R));

		S_HR = W_R(1) + a_R;
		S_TR = S_STAR + a_star_R;
	}
}





Eigen::Vector3d exact_rs_stiffenedgas :: sample_solution (const Eigen::Vector3d& W_L, const Eigen::Vector3d& W_R, double S)
{
	Eigen::Vector3d W (3);

	
	// Find appropriate part of solution and return primitives

	if (S < S_STAR)
	{
		// To the left of the contact

		if (P_STAR > W_L(2))
		{
			// Left shock
			
			if (S < S_L)
			{
				W = W_L;
			}
			else
			{
				W(0) = rho_star_L;
				W(1) = S_STAR;
				W(2) = P_STAR;
			}
		}
		else
		{
			// Left rarefaction
			
			if (S < S_HL)
			{
				W = W_L;
			}
			else
			{
				if (S > S_TL)
				{
					W(0) = rho_star_L;
					W(1) = S_STAR;
					W(2) = P_STAR;
				}
				else
				{
					set_left_rarefaction_fan_state(W_L, S, W);
				}
			}
		}
	}
	else
	{
		// To the right of the contact

		if (P_STAR > W_R(2))
		{
			// Right shock
			
			if (S > S_R)
			{
				W = W_R;
			}
			else
			{
				W(0) = rho_star_R;
				W(1) = S_STAR;
				W(2) = P_STAR;
			}
		}
		else
		{
			// Right rarefaction
			
			if (S > S_HR)
			{
				W = W_R;
			}
			else
			{
				if (S < S_TR)
				{
					W(0) = rho_star_R;
					W(1) = S_STAR;
					W(2) = P_STAR;
				}
				else
				{
					set_right_rarefaction_fan_state(W_R, S, W);
				}
			}
		}
	}

	return W;
}





double exact_rs_stiffenedgas :: find_p_star_newtonraphson (
	
	const double rho_L,
	const double u_L,
	const double p_L,
	const double rho_R,
	const double u_R,
	const double p_R
)
{
	double p_star;
	double p_star_next;
	const double TOL = 1e-6;


	// First we set the initial guess for p_star using a simple mean-value approximation
		
	p_star_next = 0.5*(p_L+p_R);
	int n = 0;
	
	
	// Now use the Newton-Raphson algorithm

	do {
		p_star = p_star_next;

		p_star_next = p_star - total_pressure_function(p_star,rho_L,u_L,p_L,rho_R,u_R,p_R)/total_pressure_function_deriv(p_star,rho_L,p_L,rho_R,p_R);
		
		p_star_next = std::max(p_star_next, TOL);
		
		n++;

	} while ((fabs(p_star_next - p_star)/(0.5*(p_star+p_star_next)) > TOL) && n < 100);
	
	if (n == 100) p_star_next = 0.5*(p_L+p_R);

	return p_star_next;
}






double exact_rs_stiffenedgas :: total_pressure_function (

	const double p_star,
	const double rho_L,
	const double u_L,
	const double p_L,
	const double rho_R,
	const double u_R,
	const double p_R
)
{
	return	f(p_star, rho_L, p_L, gamma_L, pinf_L)
		+ f(p_star, rho_R, p_R, gamma_R, pinf_R)
		+ u_R - u_L;
}




double exact_rs_stiffenedgas :: total_pressure_function_deriv (

	const double p_star,
	const double rho_L,
	const double p_L,
	const double rho_R,
	const double p_R
)
{
	return 	f_deriv (p_star, rho_L, p_L, gamma_L, pinf_L)
		+ f_deriv (p_star, rho_R, p_R, gamma_R, pinf_R);
}






double exact_rs_stiffenedgas :: f (double p_star, double rho, double p, double gamma, double pinf)
{
	if (p_star > p)
	{
		return (p_star - p)/Q_K(p_star, rho, p, gamma, pinf);
	}
	else
	{
		return (2.0*a(rho,p,gamma,pinf)/(gamma-1.0))*(std::pow((p_star + pinf)/(p + pinf), (gamma-1.0)/(2.0*gamma)) - 1.0);
	}
}




double exact_rs_stiffenedgas :: f_deriv (double p_star, double rho, double p, double gamma, double pinf)
{
	double A = 2.0/((gamma+1.0)*rho);
	double B = (p+pinf)*(gamma-1.0)/(gamma+1.0);

	if (p_star > p)
	{
		return sqrt(A/(B+p_star+pinf))*(1.0 - ((p_star-p)/(2.0*(B+p_star+pinf))));
	}
	else
	{
		return (1.0/(rho*a(rho,p,gamma,pinf)))*std::pow((p_star+pinf)/(p+pinf), -(gamma+1.0)/(2.0*gamma));
	}
}





void exact_rs_stiffenedgas :: set_left_rarefaction_fan_state (const Eigen::Vector3d& W_L, double S, Eigen::Vector3d& W)
{
	double a_L = a(W_L(0),W_L(2),gamma_L,pinf_L);
	W(0) = W_L(0)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), 2.0/(gamma_L - 1.0));
	W(1) = (2.0/(gamma_L+1.0))*(a_L + S + ((gamma_L-1.0)/2.0)*W_L(1));
	W(2) = (W_L(2) + pinf_L)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), (2.0*gamma_L)/(gamma_L-1.0)) - pinf_L;
}





void exact_rs_stiffenedgas :: set_right_rarefaction_fan_state (const Eigen::Vector3d& W_R, double S, Eigen::Vector3d& W)
{
	double a_R = a(W_R(0),W_R(2),gamma_R,pinf_R);
	W(0) = W_R(0)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), 2.0/(gamma_R - 1.0));
	W(1) = (2.0/(gamma_R+1.0))*(- a_R + S + ((gamma_R-1.0)/2.0)*W_R(1));
	W(2) = (W_R(2) + pinf_R)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), (2.0*gamma_R)/(gamma_R-1.0)) - pinf_R;
}







double exact_rs_stiffenedgas :: Q_K (double p_star, double rho, double p, double gamma, double pinf)
{
	double A = 2.0/((gamma+1.0)*rho);
	double B = (p+pinf)*(gamma-1.0)/(gamma+1.0);
	return sqrt((p_star+pinf+B)/A);
}
	



double exact_rs_stiffenedgas :: a (double rho, double p, double gamma, double pinf)
{
	return sqrt(gamma*((p+pinf)/rho));
}





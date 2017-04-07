/*
 *	This file contains a class to solve for the exact solution of the Riemann Problem for the one dimensional Euler
 *	equations with stiffened gas equation of state
 *
 *	Author: Murray Cutforth
 *	Date:	07/04/2017
 */



#ifndef EXACT_RS_STIFFGAS_H
#define EXACT_RS_STIFFGAS_H


#include <blitz/array.h>




class exact_rs_stiffenedgas {

	public:
	
	const double gamma_L;
	const double gamma_R;
	const double pinf_L;
	const double pinf_R;
	
	double S_STAR;
	double P_STAR;
	double rho_star_L;
	double rho_star_R;
	
	double S_L;
	double S_R;
	double S_HL;
	double S_TL;
	double S_HR;
	double S_TR;

	exact_rs_stiffenedgas (double gamma_L, double gamma_R, double pinf_L, double pinf_R);
	
	
	// Functions used to generate exact solutions to Riemann problems

	void solve_RP (const blitz::Array<double,1>& W_L, const blitz::Array<double,1>& W_R);

	blitz::Array<double,1> sample_solution (const blitz::Array<double,1>& W_L, const blitz::Array<double,1>& W_R, double S);
	
	
	// Function called from 2D Godunov-type methods
	
	void celledge_primitives_2D (
		
		const blitz::Array<double,1>& W_L, 
		const blitz::Array<double,1>& W_R,
		blitz::Array<double,1>& soln
	);


	
	// Functions used to solve for p_star iteratively

	double find_p_star_newtonraphson (
	
		const double rho_L,
		const double u_L,
		const double p_L,
		const double rho_R,
		const double u_R,
		const double p_R
	);

	double total_pressure_function (

		const double p_star,
		const double rho_L,
		const double u_L,
		const double p_L,
		const double rho_R,
		const double u_R,
		const double p_R
	);

	double total_pressure_function_deriv (

		const double p_star,
		const double rho_L,
		const double p_L,
		const double rho_R,
		const double p_R
	);

	double f (double p_star, double rho, double p, double gamma, double pinf);

	double f_deriv (double p_star, double rho, double p, double gamma, double pinf);



	// Functions to find the state inside a rarefaction fan

	void set_left_rarefaction_fan_state (const blitz::Array<double,1>& W_L, double S, blitz::Array<double,1>& W);

	void set_right_rarefaction_fan_state (const blitz::Array<double,1>& W_R, double S, blitz::Array<double,1>& W);



	// Misc functions

	double Q_K (double p_star, double rho, double p, double gamma, double pinf);

	

	// Equation of state functions

	double a (double rho, double p, double gamma, double pinf);

};



#endif

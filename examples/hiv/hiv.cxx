//////////////////////////////////////////////////////////////////////////
////////////////              hiv.cxx                /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Bryson maximum range problem      ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:     Bryson and Ho (1975)              ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

typedef struct{
    double s;
    double mu1;
    double mu2;
    double mu3;
    double r;
    double Tmax;
    double k;
    double N;
    double B;
} Mydata;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    Mydata* md = (Mydata*) workspace->problem->user_data;

    double B = md->B;

    adouble T  = states[ 0 ];
    adouble Ti = states[ 1 ];
    adouble V  = states[ 2 ];
    adouble u  = controls[ 0 ];

    adouble result = -(T-(0.5*B*(1.0-u)*(1.0-u)));

    return result;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble Tdot, Tidot, Vdot;

   adouble T = states[ 0 ];
   adouble Ti = states[ 1 ];
   adouble V = states[ 2 ];

   adouble u = controls[ 0 ];

   Mydata* md = (Mydata*) workspace->problem->user_data;
   double s = md->s;
   double mu1 = md->mu1;
   double mu2 = md->mu2;
   double mu3 = md->mu3;
   double r = md->r;
   double Tmax = md->Tmax;
   double k = md->k;
   double N = md->N;

   Tdot = (s/(1.0+V)-mu1*T+r*T*(1.0-(T+Ti)/(Tmax)))-u*k*V*T;
   Tidot = u*k*V*T-mu2*Ti;
   Vdot = N*mu2*Ti-mu3*T*V;

   derivatives[ 0 ] = Tdot;
   derivatives[ 1 ] = Tidot;
   derivatives[ 2 ] = Vdot;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble T0 = initial_states[ 0 ];
   adouble Ti0 = initial_states[ 1 ];
   adouble V0 = initial_states[ 2 ];

   e[ 0 ] = T0;
   e[ 1 ] = Ti0;
   e[ 2 ] = V0;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    Mydata* md = (Mydata*) malloc(sizeof(Mydata));
    problem.user_data = (void*) md;
    double s = 10.0;
    double mu1 = 0.02;
    double mu2 = 0.5;
    double mu3 = 4.4;
    double r = 0.03;
    double Tmax = 1500.;
    double k = 0.000024;
    double N = 300.;
    double B = 10.0;

    md->s = s;
    md->mu1 = mu1;
    md->mu2 = mu2;
    md->mu3 = mu3;
    md->r = r;
    md->Tmax = Tmax;
    md->k = k;
    md->N = N;
    md->B = B;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "HIV Treatment Problem";
    problem.outfilename         = "hiv.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages           = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 3;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes       << 100;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double t0=0.0;
    double tf=100.0;

    double TL = 0.0;
    double TiL = 0.0;
    double VL = 0.0;
    double TU = Tmax;
    double TiU = PSOPT::inf;
    double VU = PSOPT::inf;

    double uL = 0.0;
    double uU = 1.0;

    double T0 = 806.4;
    double Ti0 = 0.04;
    double V0 = 1.5;

    problem.phases(1).bounds.lower.states(0) = TL;
    problem.phases(1).bounds.lower.states(1) = TiL;
    problem.phases(1).bounds.lower.states(2) = VL;

    problem.phases(1).bounds.upper.states(0) = TU;
    problem.phases(1).bounds.upper.states(1) = TiU;
    problem.phases(1).bounds.upper.states(2) = VU;


    problem.phases(1).bounds.lower.controls(0) = uL;
    problem.phases(1).bounds.upper.controls(0) = uU;


    problem.phases(1).bounds.lower.events(0) = T0;
    problem.phases(1).bounds.lower.events(1) = Ti0;
    problem.phases(1).bounds.lower.events(2) = V0;

    problem.phases(1).bounds.upper.events(0) = T0;
    problem.phases(1).bounds.upper.events(1) = Ti0;
    problem.phases(1).bounds.upper.events(2) = V0;


    problem.phases(1).bounds.lower.StartTime    = t0;
    problem.phases(1).bounds.upper.StartTime    = t0;

    problem.phases(1).bounds.lower.EndTime      = tf;
    problem.phases(1).bounds.upper.EndTime      = tf;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             = &dae;
    problem.events 		    = &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    = ones(nstates,nnodes);

    x_guess.row(0)  = T0*ones(1,nnodes);
    x_guess.row(1)  = Ti0*ones(1,nnodes);
    x_guess.row(2)  = V0*ones(1,nnodes);

    problem.phases(1).guess.controls       = 0.5*ones(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(t0,tf,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.mesh_refinement             = "automatic";
//    algorithm.collocation_method          = "trapezoidal";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"lambda.dat");
    Save(H,"H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name+": states", "time (s)", "states","T Ti V");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u");

    plot(t,x,problem.name+": states", "time (s)", "states","T Ti V",
                             "pdf", "hiv_states.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u",
                             "pdf", "hiv_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

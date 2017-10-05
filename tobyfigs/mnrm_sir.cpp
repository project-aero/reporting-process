/*ABOUT: This is an implementation of the MNRM as found in Anderson (2008).
*/

//headers and declarations
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <gsl/gsl_rng.h> 
#include <list>

using namespace std;


/**** System Parameters ****/

// System size:
double N = 10000; 

// Number of reaction channels:
const int a_no = 5; 
// Number of species:
const int s_no = 2;
//Number of output channels to datafile:
const int no_output = 6;
// Initial relaxation time:
const double temin = 10;
// Output time-step size:
const double dte= 1.0; //.0/2048.0;



/****  Model Parameters ****/
//Units: weeks
double R0 = 1.0; //final value of R0.
double gamm = 1.0;
const double zeta = 1.0;
const double mu = 1./(70.*52.);
const double nu = 1./(70.*52.);
double beta = R0*(gamm +mu);
const double p0 = 0.0;
const double tvacc = temin;



/**** Function definitions ****/

// Reaction rates:
void reactions_update(double  n[s_no],  double a[a_no], double vaccine_uptake)
				{			
				//Infection:	
				a[0] = (n[0]/N)*(beta*vaccine_uptake*n[1] + zeta);
				//Recovery:
				a[1] = gamm*n[1] ;
				//Birth:
				a[2] = nu*N;
				//Death S:
				a[3] = mu*n[0];
				//Death I:
				a[4] = mu*n[1];
				}
// Stoichiometric Matrix:
const int v[s_no*a_no] = {-1,1,
						   0,-1,
						   1,0,
						   -1,0,
						   0,-1};

//Vaccine uptake function:
double vaccine_uptake( double t, double initial_vaccine, double rate_of_vaccination){return( initial_vaccine + rate_of_vaccination*t);}
//Deterministic fixed point:
void fixed_point(double nstar[s_no], double uptake);
//Reset initial conditions between runs:
void reset_initial_conditions( double n[s_no], double uptake);




ofstream outfreqrunning;

int main(int argc, const char* argv[])
{



	double runs = 1;
	int L = 52*100;
	double null_param = 0;
	srand (time(NULL));
	
	
	if(argc > 1)runs = atoi(argv[1]);
	else cerr << "# No value for runs given, default used: runs = 1" << endl;
	if(argc > 2)N = atoi(argv[2]);
	else cerr << "# No value for N given, default used: N = 1000" << endl;
	if(argc > 3)L = 52*atoi(argv[3]); // Number of output time-steps
	else cerr << "# No value for L given, default used: L = 100" << endl;
	if(argc > 4) srand (atoi(argv[4])); // seed for rng
	else cerr << "# No seed for rng given, default used: srand (time(NULL))" << endl;
	if(argc > 5){if(strncmp(argv[5],"-n",2) == 0) null_param = 1;} // seed for rng
	else cerr << "# No value for scenario, default used: -e (emerging)" << endl;
	if(argc > 6){if(strncmp(argv[6],"-monthly",2) == 0){gamm= 7./30.;}} // seed for rng
	else cerr << "# No value for infectious period, default used: -weekly" << endl;

	beta = R0*(gamm +mu);
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);	
	gsl_rng_set (rng, rand());

		
	// Final time:
	double temax = temin + (double)L*dte;
	double rate_of_uptake = 1/(dte*(double)L);

	const double rtN = sqrt(N);
	if(argc > 4) srand( atoi(argv[4]));
	else srand ( time(NULL) );
	// Output filenames:
	char model_filename [50], ps_filename [50], ac_filename[50], model_type [50];
	sprintf (model_type, "sis");
	sprintf (model_filename, "%s", model_type);




	//Deterministic fixed point:
	double nstar[s_no] = {0, 0};


/**** System variables: ****/
	double t = 0;
	double a[a_no];
	int q = 0;
	int nu;
	double auto_y = 0;

	double dt;
	double te = temin;
	double ps_m[L], var_m[L], xi_m[L], ns[L],  comp[2*L]; 

	// Internal Poisson processes, internal clocks, next-firing times:
	double P[a_no], T[a_no], D[a_no];

	int r = 0;
	//Level of vaccine uptake at time 0:
	double vaccine_uptake = p0;
	fixed_point(nstar, vaccine_uptake);
	//cout << "# fixed point: " << nstar[0] << endl;	
	double n[s_no] = {N-1, 1};
	double output[L][no_output];

while(r<runs)
	{
	//Reset initial conditions:
	reset_initial_conditions(n, vaccine_uptake);
	if(null_param == 1) vaccine_uptake = 0.0;
	else vaccine_uptake = 0.0;
	t = 0;
	te = temin;
	q = 0;
	int cases = 0;


	

	//Initialisation of internal Poisson processes and internal clocks:
	for(int k = 0; k < a_no; k++){ P[k] = -log(gsl_rng_uniform_pos(rng)); T[k] = 0;}
	
	//Array of fluctuations about the deterministic fixed point:
	while(q < L)
		{

		//cout << t << "\t" << n[0] << "\t" << n[1] << "\t" << vaccine_uptake << endl;
		if(t>=te){
			while (te <= t && q < L){
				fixed_point(nstar, vaccine_uptake);
				xi_m[q] = (n[1] - nstar[1])/rtN;
				output[q][0] = t; //q*dte;
				output[q][1] = n[0];
				output[q][2] = n[1];
				output[q][3] += n[1]/runs ; //
				output[q][4] += n[1]*n[1]/runs;
				output[q][5] = beta*vaccine_uptake;
				cout << te-temin << " " << n[1] << " " << beta*vaccine_uptake/(gamm +mu) << " " << N << " " << L << " " << r << endl; //" " << cases << endl;
				cases = 0;
				//for( int i =0; i < no_output; i++)cout << output[q][0] << "\t";
				//cout << endl;
				q++;
				te += dte;

				}
			}

	// Updating the reaction rates:
		reactions_update(n,  a,  vaccine_uptake);
//&& a[k] > 10e-20
	// Calculating which reaction fires next:		
		dt = temax;
		for(int k = 0; k < a_no; k++){D[k] = (P[k] - T[k])/a[k]; if(D[k] <= dt && a[k] > 10e-20){dt = D[k]; nu = k;}}


	// Updating the internal Possion process and system state according to the reaction with fired:
		t = t+ dt;
		if( nu < a_no){ 
			P[nu] -= log(gsl_rng_uniform_pos (rng) );
			for(int i = 0; i < s_no; i++) n[i] += v[i+s_no*nu];}
		if(nu == 1) cases++;

	// Internal clocks are updated:
		for(int k = 0; k < a_no; k++)T[k] += dt*a[k];		

		
	 //If t < tvacc, decrease p:
		if(null_param != 1 && t > tvacc && vaccine_uptake + dt*rate_of_uptake  < 1.0) {vaccine_uptake += dt*rate_of_uptake;}

		}

	// Condition for ending a run and operations performed at the end of a run:
	if(q==L){		

		r++;	
		q = 0;
	//	cout << "# " << r << endl;

			sprintf (ps_filename, "./data/ts/ts_r_%d.txt",r);
			outfreqrunning.open(ps_filename);
			outfreqrunning << "# Parameters: L = " << L << " dte = " << dte << " N = " << N << " ys = " << nstar[1]/N << endl;
			for (int i = 0; i < L; i++){for (int j = 0; j < no_output; j++){    outfreqrunning << output[i][j] << "\t"; } outfreqrunning << endl;} 
			outfreqrunning.close();




		}
	}


return(0);
}		


	
////////////////////////////////////////////////////////////////////////////////

void reset_initial_conditions( double n[s_no], double uptake)
		{
		uptake = 0.0; 
		n[0] = N -1;
		n[1] = 1;
		}


void fixed_point(double nstar[s_no], double uptake)
		{
		double RootPart = sqrt( ( uptake*R0*mu - R0*mu - mu-zeta/N)* ( uptake*R0*mu - R0*mu - mu-zeta/N)  -4*R0*mu*(mu-uptake*mu) ) ;
		nstar[0] = (N/(2*R0*mu)) *  ( zeta/N + mu + R0*mu -uptake* R0 * mu - RootPart) ;
		nstar[1] = (N/(2*gamm + 2*mu))*(   mu*(1- uptake) - (zeta/N+ mu)/R0 + (1/R0)*RootPart);

		//nstar[1] = N*mu*((1-uptake)*R0 -1 )/beta;
		}



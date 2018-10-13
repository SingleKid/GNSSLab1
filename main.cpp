#pragma warning(disable : 4996) // for VC using standard i-o

#define _USE_MATH_DEFINES       // for using M_PI M_PI_2 defined in math.h
#include <math.h>
#include <stdio.h>


#define OBS_FILE "RemoteL1L2.obs"  // file name of remote station
#define SAT_FILE "Satellites.sat"  // file name of satellite info
#define STA_FILE "BaseL1L2.obs"    // file name of base station

#include "MatC.h"                  // Matrix lib

#define GPS_SAT_NUM 32             // number of GPS satellite
#define CHANNEL_NUM 12             // number of channel
#define LS_MAX_ITER 20             // max number of iteration in Least Square
#define LS_CONV_THRES 0.00001      // the threshold for convergence
#define SEC_OF_WEEK 604800         // seconds of one week

#define OBS_SIG0 1                 // standard deviation for pseudorange observation from Zenith (in meters)
#define ELEV_INVALID -1            // invalid value for elevation
#define ELEV_THRES 1              // elevation threshold

double STA_COOR[3] { -1625352.17084393, -3653483.75114927, 4953733.86925805 }; // accurate coordinates of base station

// For the processing modes: non-difference, single-difference, double-difference.
enum MODE {
	SPP,
	DSPP,
	DDSPP,

};

MODE now = DDSPP;    // mode for now.
int sat_limit = 5; // limitation of using satellites in one epoch

// a frame of data from obs file
struct obs_epoch {
	double PRN;    // satellite number
	double Time;   // observation time
	double C1;     // L1 pseudorange
	double L1;     // L1 carrier phase
	double D1;     // L1 doppler
	double L2;     // L2 carrier phase
};

// a frame of data from sat file
struct sat_epoch {
	double PRN;    // satellite number
	double Time;   // observation time
	double P[3];   // satellite position
	double V[3];   // satellite velocity
};

FILE * ofp = fopen(OBS_FILE, "rb");  // obs file handle
FILE * nfp = fopen(SAT_FILE, "rb");  // sat file handle

FILE * sfp = (now == DSPP || now == DDSPP )? fopen(STA_FILE, "rb") : NULL; // station file handle

obs_epoch OBS[CHANNEL_NUM];  // buffer for obs input data
sat_epoch SAT[CHANNEL_NUM];  // buffer for sat input data
obs_epoch STA[CHANNEL_NUM];  // buffer for obs station input data

bool obs_available[CHANNEL_NUM]; // flag for available obs
bool sat_available[CHANNEL_NUM]; // flag for available sat
bool sta_available[CHANNEL_NUM]; // flag for available sta

bool solve_available[CHANNEL_NUM]; // flag for satellites available for solving


obs_epoch * S_OBS[CHANNEL_NUM];    // buffer of obs data while solving
sat_epoch * S_SAT[CHANNEL_NUM];    // buffer of sat data while solving
obs_epoch * S_STA[CHANNEL_NUM];    // buffer of station data while solving

int sat_num = 0;                   // number of satellites in this epoch
double current_time;               // current gps time in this epoch

double solution[4]{ 0,0,0,0 };     // solution buffer, X, Y, Z, cdT
                                   // please note solution[3] is cdT, not dT. So A(i,4) = 1, not c.

bool solution_available = false;   // solution avaiablility flag

double sat_elevation[GPS_SAT_NUM]; // elevation result of satellites


// This function is used to check the input obs data from files/streams.
// it refuses observations with any parameters invalid or time not correct.
bool obs_check(obs_epoch * OBS)
{
	current_time = OBS->Time;

	// for every channel in buffer
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		try {
			obs_epoch * ref = OBS + i;
			// checking its prn
			if (ref->PRN <= 0 || ref->PRN >= GPS_SAT_NUM)
				throw 1;
			// checking its time
			else if (ref->Time < 0 || ref->Time >= SEC_OF_WEEK)
				throw 1;
			// checking its basic observation
			else if (ref->C1 <= 0)
				throw 1;
			// checking it's time again, for synchronizing
			else if(ref->Time != current_time)
				throw 1;

			obs_available[i] = true;
			
		}
		catch (...)
		{
			// set the flag to avoid using it.
			obs_available[i] = false;
		}
	}
}

// This function is to calculate the distance between two vectors (or coorindate locations)
inline double _fastcall distance(double * p1, double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

// This function is to calculate latitude, longitude, altitude from XYZ in WGS84
void XYZ2BLH(double * XYZ, double * BLH)
{
	const static double a = 6378137.0;
	const static double F = 1.0 / 298.257223563;
	double e2, Z, dZ, ZdZ, r, sinb, N, x2y2;
	int iter;
	iter = 0;
	r = 0.0;
	N = 0.0;
	sinb = 0.0;
	e2 = 2 * F - F * F;
	x2y2 = XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1];
	dZ = e2 * XYZ[2];
	do
	{
		Z = dZ;
		ZdZ = Z + XYZ[2];
		r = x2y2 + ZdZ*ZdZ;
		sinb = ZdZ / sqrt(r);
		N = a / sqrt(1 - e2*sinb*sinb);
		dZ = N * e2 * sinb;
		iter = iter + 1;
	} while ((iter <= 10) && (fabs(dZ - Z) > 1E-8));
	BLH[0] = atan2(XYZ[1], XYZ[0]);
	BLH[1] = atan2(ZdZ, sqrt(x2y2));
	BLH[2] = sqrt(x2y2 + ZdZ*ZdZ) - N;
}

// To get the rotation matrix from XYZ to ENU, using the latitude and longitude of the local station.
void get_matrix_T(Matrix *& E, double * BLH)
{
	E = malloc_mat(3, 3);

	double sinp = sin(BLH[0]), cosp = cos(BLH[0]), sinl = sin(BLH[1]), cosl = cos(BLH[1]);
	E->data[0][0] = -sinl;       E->data[0][1] = cosl;        E->data[0][2] = 0.0;
	E->data[1][0] = -sinp*cosl;  E->data[1][1] = -sinp*sinl;  E->data[1][2] = cosp;
	E->data[2][0] = cosp*cosl;   E->data[2][1] = cosp*sinl;   E->data[2][2] = sinp;

}

// To get the satellite elevation from the positions of satellites and user.
double elev(double * sat_pos, double * user_pos)
{

	double dpos[3] = { 0 };
	double ori[3]{ 0,0,0 };
	dpos[0] = sat_pos[0] - user_pos[0];
	dpos[1] = sat_pos[1] - user_pos[1];
	dpos[2] = sat_pos[2] - user_pos[2];

	double user_distance_to_earth = distance(user_pos, ori);

	double mod = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
	if (fabs(user_distance_to_earth * mod < 1.0)) {
		return M_PI_2;
	}
	else {
		double m = dpos[0] * user_pos[0] + dpos[1] * user_pos[1] + dpos[2] * user_pos[2];
		double n = m / (mod * user_distance_to_earth);
		return M_PI_2 - acos(n);
	}
}

// Check the input data from satellite file, avoid using invalid satellites.
bool sat_check()
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		try {
			sat_epoch * ref = SAT + i;
			int PRN = ref->PRN;
			
			// for prn
			if (ref->PRN <= 0 || ref->PRN >= GPS_SAT_NUM)
				throw 1;
			// for time
			else if (ref->Time < 0 || ref->Time >= SEC_OF_WEEK)
				throw 1;
			//for synchronizing
			else if (ref->Time != current_time)
				throw 1;

			sat_elevation[PRN] = ELEV_INVALID;

			// if the reciever solution is available, then get the elevation
			if (solution_available) {
				double e = elev(ref->P, solution);
				// if elevation is greater than the threshold, we can use it.
				if (e >= ELEV_THRES * M_PI / 180.0) {
					sat_available[i] = true;
					sat_elevation[PRN] = e;
				}
				// otherwise we cannot.
				else
					sat_available[i] = false;
			}
			// if no solution is available, then use every satallite.
			else {
				sat_available[i] = true;
			}
		}
		catch (...)
		{
			sat_available[i] = false;
		}
	}
}

// check all the parameters again and combine them into use, build the second order buffers: 
// use OBS, SAT, STA
// use obs_available, sat_available,sta_available
// build S_OBS, S_SAT, S_STA
// build solve_available, sat_num
void overall_check()
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		if (!obs_available[i]) continue; // if obs not available, this channel is done.

		int prn = OBS[i].PRN;            // get the prn number for multiple uses inside this function, to make program faster
		for (int j = 0; j < CHANNEL_NUM; j++)  // search the prn in the satellite list
		{
			if (SAT[j].PRN == prn && sat_available[j]) // if found
			{
				if (now == SPP) {  // if un-differenced mode, there is no need of station observation, parse the data and done.
					S_OBS[sat_num] = OBS + i; 
					S_SAT[sat_num] = SAT + j;
					S_STA[sat_num] = NULL; // no need of STA 
					sat_num++;
					solve_available[i] = true;
				}
				else if (now == DSPP || now == DDSPP) {
					for (int k = 0; k < CHANNEL_NUM; k++) // try to find prn in station observations.
					{
						if (STA[k].PRN == prn && sta_available[k])
						{
							S_OBS[sat_num] = OBS + i;
							S_SAT[sat_num] = SAT + j;
							S_STA[sat_num] = STA + k;
							sat_num++;
							solve_available[i] = true;
							break;
						}
					}
				}
				break;
			}
			else {
				solve_available[i] = false; // if fail, this prn is not solve available.
			}
		}
	}
}

Matrix * residuals; // get the residuals from solver to get it extracted in file.
double DOPs[5]{ 0,0,0,0,0 }; // ENVHP, DOPS.

void overall_reset()  // do this function after each epch to get everything re-initialized.
{
	sat_num = 0;
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		sat_available[i] = false;
		obs_available[i] = false;
		current_time = 0;
	}
	free_mat(residuals);  // clean the buffer of residuals.
}

int DD_ref_sat = 0;   // this is the reference prn for double-differencing mode.

// This function is to apply Least Square method in SPP.
bool solve()
{
	// epoch initialization
	Matrix * L  = malloc_mat(sat_num,       1);       // ***OBERVATION  VECTOR***    which is also used as ***MISCLOSURE VECTOR***
	Matrix * A  = malloc_mat(sat_num,       4);       // ***DESIGN      MATRIX***
	Matrix * Cl = malloc_mat(sat_num, sat_num);       // ***COVARIANCE  MATRIX***
	Matrix * r  = malloc_mat(sat_num,       1);       // ***RESIDUAL    VECTOR***
	Matrix * Q  = malloc_mat(      4,       4);       // ***COFACTOR    MATRIX***
	Matrix * δ  = malloc_mat(      4,       1);       // ***UNKNOWN     VECTOR***

	// for double-differencing mode
	Matrix * D = NULL;                                // *** HOW TO CALL THIS ***
	if (now == DDSPP)
	{
		// try find the reference satellite which is just the highest one.
		// But if the solution is unavailable, the elevations can not be obtained.
		// which means we just use the first one in channels.
		double ref_ele = sat_elevation[(int)(S_SAT[0]->PRN)]; 
		DD_ref_sat = 0;

		if (solution_available) { // if we have located
			// find the highest.
			for (int i = 1; i < sat_num; i++)
			{
				if (ref_ele < sat_elevation[(int)(S_SAT[i]->PRN)])
				{
					ref_ele = sat_elevation[(int)(S_SAT[i]->PRN)];
					DD_ref_sat = i;
				}
			}
		}

		// Build the D matrix using voted reference satellite.
		D = malloc_mat(sat_num - 1, sat_num);
		for (int i = 0; i < sat_num - 1; i++) 
		{
			D->data[i][i + (i >= DD_ref_sat ? 1 : 0)] = 1;
			D->data[i][DD_ref_sat] = -1;
		}
	}

	// some buffers to build the A matrix and L matrix
	double * DX0 = (double*)alloca(sat_num * sizeof(double));  // Xs - Xr
	double * DY0 = (double*)alloca(sat_num * sizeof(double));  // Ys - Yr
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));  // Zs - Zr
	double * S   = (double*)alloca(sat_num * sizeof(double));  // Approx distance between s and r : ρ

	// for between reciever single differencing, another buffer
	double * S2 = (double*)alloca(sat_num * sizeof(double));   // Approx distance between s and station :ρ2

	// for iteration
	double last_solution[4] = { 0,0,0,0 }; //last solution


	for (int i = 0; i < LS_MAX_ITER; i++) // do iteration with a max iteration number
	{
		memcpy(last_solution, solution, sizeof(double) * 4); // copy the solution

		// before getting A and L, buffers must be filled.
		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_SAT[j]->P[0] - solution[0];
			DY0[j] = S_SAT[j]->P[1] - solution[1];
			DZ0[j] = S_SAT[j]->P[2] - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);

			// for between reciever single differencing
			if(now == DSPP || now == DDSPP)
				S2[j] = distance(STA_COOR, S_SAT[j]->P);
		}

		// get A, L, Cl matrices.
		for (int j = 0; j < sat_num; j++)
		{
			if (sat_elevation[(int)(S_SAT[j]->PRN)] == ELEV_INVALID) // if solution not avaiable, set the covariance matrix using constant
				Cl->data[j][j] = OBS_SIG0 * OBS_SIG0 / sin(ELEV_THRES * M_PI / 180.0);
			else                                                     // if solution available, set the covariance matrix using elevations.
				Cl->data[j][j] =  OBS_SIG0 * OBS_SIG0 / sin(sat_elevation[(int)(S_SAT[j]->PRN)]);

			// for between reciever single differencing
			if (now == DSPP || now == DDSPP)
				L->data[j][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j] - S[j] - solution[3];

			// for common spp
			else if (now == SPP) L->data[j][0] = S_OBS[j]->C1 - S[j] - solution[3];

			// A matrix
			A->data[j][0] = -DX0[j] / S[j];  // for X
			A->data[j][1] = -DY0[j] / S[j];  // for Y
			A->data[j][2] = -DZ0[j] / S[j];  // for Z
			A->data[j][3] = 1;               // for cdt
		}


		// double differencing, get the new L, A, Cl, using D matrix.
		if (now == DDSPP)
		{
			Matrix * La = L;
			Matrix * Aa = A;
			Matrix * Cla = Cl, * temp = NULL, *Dt = NULL;

			A = NULL;
			L = NULL;
			Cl = NULL;
			mat_multiply(D, La, L); // Ld = D * L
			mat_multiply(D, Aa, A); // Ad = D * A

			mat_multiply(D, Cla, temp); // Cld = D * Cl * D'
			mat_trans(D, Dt);
			mat_multiply(temp, Dt, Cl);

			A->cols = 3;            // Cancel the last colomn in A (which is the clock offset). 
			                        // Please note this code is not memory-safe for my MatC lib but a temporary solution.

			// clean memory
			free_mat(Aa);
			free_mat(La);
			free_mat(Cla);
			free_mat(temp);
			free_mat(Dt);
		}

		LMS(L, A, Cl, δ, Q, r);   // Do LS!!!!!!!!
		                           // Please note the output Q = inv(A'PA), but we need inv(A'A) for DOPS.

		if (now == DDSPP) {
			for (int j = 0; j < 3; j++)
				solution[j] += δ->data[j][0];
			solution[3] = 0; // clock offset is zero in this mode.
		}
		else
			for (int j = 0; j < 4; j++)
				solution[j] += δ->data[j][0];

		// If converged
		if (distance(last_solution, solution, 3) <= LS_CONV_THRES) {

			// job done
			free_mat(Q);
			Matrix * T = NULL, * Tt = NULL, * temp1 = NULL, * Qn = NULL;

			// get Q = inv(A'A)
			Matrix * At = NULL, *AtA = NULL, *Q = NULL;
			mat_trans(A, At);
			mat_multiply(At, A, AtA);
			mat_inv(AtA, Q);
			double blh[3];
			residuals = r;

			double x2 = Q->data[0][0] * Q->data[0][0];
			double y2 = Q->data[1][1] * Q->data[1][1];
			double z2 = Q->data[2][2] * Q->data[2][2];
			//double t2 = Q->data[3][3] * Q->data[3][3];

			//ENVHP DOPS calculation
			Q->cols = 3;
			Q->rows = 3;
			DOPs[4] = sqrt(x2 + y2 + z2);// PDOP
			XYZ2BLH(solution, blh);
			get_matrix_T(T, blh);

			// get Qx = TQT'
			mat_inv(T, Tt);
			mat_multiply(T, Q, temp1);
			mat_multiply(temp1, Tt, Qn);

			double e2 = Qn->data[0][0] * Qn->data[0][0];
			double n2 = Qn->data[1][1] * Qn->data[1][1];
			double u2 = Qn->data[2][2] * Qn->data[2][2];

			DOPs[0] = sqrt(e2);        // EDOP
			DOPs[1] = sqrt(n2);        // NDOP
			DOPs[2] = sqrt(u2);        // VDOP
			DOPs[3] = sqrt(e2 + n2);   // HDOP

			// clean memory
			free_mat(T);
			free_mat(Tt);
			free_mat(temp1);
			free_mat(Qn);
			free_mat(At);
			free_mat(AtA);

			return true;
		}

		// if haven't converge yet, this is to reset matrix A L Cl for double differencing mode.
		if (now == DDSPP)
		{
			free_mat(A);
			free_mat(L);
			free_mat(Cl);
			L = malloc_mat(sat_num, 1);
			A = malloc_mat(sat_num, 4);
			Cl = malloc_mat(sat_num, sat_num);
		}
	}

	return false;
}


// Open all the output files.
FILE * out1    = fopen("out1.txt"  , "w");     // solution  file
FILE * logs    = fopen("log.txt"   , "w");     // residual  file
FILE * dopf    = fopen("dops.txt"  , "w");     // dops      file
FILE * satn    = fopen("satn.txt"  , "w");     // sat num   file
FILE * res2    = fopen("res2.txt"  , "w");     // elevation file
FILE * ref_sat = fopen("refsat.txt", "w");     // reference satellite for every epoch (for double-differ mode of course)

// this function is to output epoch information in files.
void output()
{
	fprintf(out1, "%lf %lf %lf %lf\n", 
		solution[0], solution[1], solution[2], solution[3]);

	if (now == DDSPP) sat_num--;
	for (int i = 0; i < GPS_SAT_NUM; i++)
	{
		for (int j = 0; j < sat_num; j++)
		{
			if (S_OBS[j]->PRN == i + 1)
			{
				fprintf(logs, "%5.3lf\t", -residuals->data[j][0]);
				fprintf(res2, "%5.3lf\t", sat_elevation[(int)(S_SAT[j]->PRN)]);
				goto next;
			}
		}
		fprintf(logs, " 0.000\t");
		fprintf(res2, " 0.000\t");
	next:;
	}

	fprintf(logs, "\n");
	fprintf(res2, "\n");


	if (now == DDSPP)
	{
		sat_num++;
		fprintf(ref_sat, "%lg\t%lf\n", S_OBS[DD_ref_sat]->PRN, sat_elevation[(int)S_OBS[DD_ref_sat]->PRN]);
	}

	fprintf(dopf, "%lf %lf %lf %lf %lf\n",
		DOPs[0], DOPs[1], DOPs[2], DOPs[3], DOPs[4]);

	fprintf(satn, "%d\n", sat_num);
}

int main()
{
	// main program 
	while (!feof(ofp) && !feof(nfp)) // execute till anyone reached EOF
	{
		if (now == DSPP || now == DDSPP) { // see wether we need station data
			if (feof(sfp))break;           // if EOF in station file, end the program
			fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);  // read station data
			obs_check(STA);                                   // check station data
		}

		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);// read observation data
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);// read satellite data

		obs_check(OBS); // check observation data
		sat_check();    // check satellite data

		overall_check(); // combine data together for processing

		// see if sat_num is larger than 4. otherwise it can't be executed.
		if (sat_num >= 4) {
			
			if (sat_num > sat_limit)sat_num = sat_limit; // this is for the limited satellite cases. a sat_limit can be set in.
			solution_available = solve();                // solve it.

			output();                                    // output it.
		}
		overall_reset(); // reset buffers and flags for another epoch.		
	}

	_fcloseall();    // save all the files.
}
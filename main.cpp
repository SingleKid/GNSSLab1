#pragma pack(1)
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>


#define OBS_FILE "obs.file"
#define SAT_FILE "sat.file"
#define STA_FILE "sta.file"

#include "MatC.h"

#define GPS_SAT_NUM 32
#define CHANNEL_NUM 12
#define LS_MAX_ITER 10
#define LS_CONV_THRES 0.0001
#define SEC_OF_WEEK 604800

#define OBS_SIG0 1 // m

enum MODE {
	SPP,
	DSPP
};
MODE now = SPP;

struct obs_epoch {
	double PRN;
	double Time;
	double C1;
	double L1;
	double D1;
	double L2;
};

struct sat_epoch {
	double PRN;
	double Time;
	double P[3];
	double V[3];
};

FILE * ofp = fopen(OBS_FILE, "rb");
FILE * nfp = fopen(SAT_FILE, "rb");

FILE * sfp = now == DSPP ? fopen(STA_FILE, "rb") : NULL;

obs_epoch OBS[CHANNEL_NUM];
sat_epoch SAT[CHANNEL_NUM];
obs_epoch STA[CHANNEL_NUM];

bool obs_available[CHANNEL_NUM];
bool sat_available[CHANNEL_NUM];
bool sta_available[CHANNEL_NUM];

bool solve_available[CHANNEL_NUM];


obs_epoch * S_OBS[CHANNEL_NUM];
sat_epoch * S_SAT[CHANNEL_NUM];
obs_epoch * S_STA[CHANNEL_NUM];

int sat_num = 0;
double current_time;

double solution[4];
bool solution_available = false;

double sat_elevation[CHANNEL_NUM];

bool obs_check(obs_epoch * OBS)
{
	current_time = OBS->Time;

	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		try {
			obs_epoch * ref = OBS + i;
			if (ref->PRN <= 0 || ref->PRN >= GPS_SAT_NUM)
				throw 1;
			else if (ref->Time < 0 || ref->Time >= SEC_OF_WEEK)
				throw 1;
			else if (ref->C1 <= 0)
				throw 1;
			else if(ref->Time != current_time)
				throw 1;

			obs_available[i] = true;
			
		}
		catch (...)
		{
			obs_available[i] = false;
		}
	}
}

inline double _fastcall distance(double * p1, double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

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

bool sat_check()
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		try {
			sat_epoch * ref = SAT + i;
			if (ref->PRN <= 0 || ref->PRN >= GPS_SAT_NUM)
				throw 1;
			else if (ref->Time < 0 || ref->Time >= SEC_OF_WEEK)
				throw 1;
			else if (ref->Time != current_time)
				throw 1;

			sat_available[i] = true;
			sat_elevation[i] = solution_available ? elev(ref->P, solution) : -1;
		}
		catch (...)
		{
			sat_available[i] = false;
		}
	}
}

void overall_check()
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		bool sta = ((now == DSPP) ? sta_available[i] : true);
		if (sat_available[i] && obs_available[i] && sta)
		{
			S_OBS[sat_num] = OBS + i;
			S_SAT[sat_num] = SAT + i;
			S_STA[sat_num] = sta ? STA + i : NULL;

			sat_num++;
			solve_available[i] = true;
		}
		else {
			solve_available[i] = false;
		}
	}
}

void overall_reset()
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		sat_available[i] = false;
		obs_available[i] = false;
		current_time = 0;
	}
}

bool solve()
{
	// epoch initialization
	Matrix * L = malloc_mat(sat_num, 1);
	Matrix * A = malloc_mat(sat_num, 4);
	Matrix * Cl = malloc_mat(sat_num, sat_num);
	Matrix * r = malloc_mat(sat_num, 1);
	Matrix * Sig = malloc_mat(4, 4);
	Matrix * δ = malloc_mat(4, 1);


	double * DX0 = (double*)alloca(sat_num);
	double * DY0 = (double*)alloca(sat_num);
	double * DZ0 = (double*)alloca(sat_num);
	double * S   = (double*)alloca(sat_num);

	double last_solution[4]; //solution


	for (int i = 0; i < LS_MAX_ITER; i++)
	{
		memcpy(last_solution, solution, sizeof(double) * 4);
		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_SAT[j]->P[0] - solution[0];
			DY0[j] = S_SAT[j]->P[1] - solution[1];
			DZ0[j] = S_SAT[j]->P[2] - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);
		}

		for (int j = 0; j < sat_num; j++)
		{
			Cl->data[j][j] = OBS_SIG0 * OBS_SIG0 / sin(sat_elevation[j]);
			L->data[j][0] = S_OBS[j]->C1 - S[j] - solution[3];
			A->data[j][0] = -DX0[j] / S[j];
			A->data[j][1] = -DY0[j] / S[j];
			A->data[j][2] = -DZ0[j] / S[j];
			A->data[j][3] = 1;
		}

		LMS(L, A, Cl, δ, Sig, r);

		for (int j = 0; j < 4; j++)
			solution[j] += δ->data[j][0];

		if (distance(last_solution, solution, 3) <= LS_CONV_THRES) {
			// job done

			break;
		}
	}
}

int main()
{

	while (!feof(ofp) && !feof(nfp))
	{
		if (now == DSPP) {
			if (feof(sfp))break;
			fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);
			obs_check(STA);
		}

		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);

		obs_check(OBS);
		sat_check();

		overall_check();
		solve();
		overall_reset();

	}
}
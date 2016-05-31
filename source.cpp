# include <iostream>
# include <fstream>
# include <string.h>
# include <algorithm>
# include <math.h>




using namespace std;






void navierstokes_semi_explicit()
{
	//inputs
	double rho = 10; double myu = 0.007; //viscosity
	double k = 100;
	double alpha = 0.8;
    int imax = 4;
	int jmax = 4;
	double L1 = 1, L2 = 1;
	double u_vel_west = 0, u_vel_east = 0, u_vel_south = 0, u_vel_north = 1;
	double v_vel_west = 0, v_vel_east = 0, v_vel_south = 0, v_vel_north = 0;
	double q_gen = 10;
	int i = 0, j = 0; double c = 300;//cP value or one in case of fluid

	//STEP-2: Geometrical Parameter and Stability criterion based time-step
	
	//double alpha = k / (rho*cp);
	double conver_crit = 0.08;
	double Dx = L1 / (imax);
	double Dy = L2 / (jmax);
	double DV = Dx*Dy;

	



	


	double ** mx1;
	mx1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mx1[j] = new double[imax];

	double ** ax1;
	ax1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ax1[j] = new double[imax];

	double ** dx1;
	dx1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dx1[j] = new double[imax];


	double ** my1;
	my1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my1[j] = new double[imax];


	double ** ay1;
	ay1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ay1[j] = new double[imax];

	double ** dy1;
	dy1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dy1[j] = new double[imax];

	double ** my1;
	my1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my1[j] = new double[imax];


	double ** mx2;
	mx2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ay1[j] = new double[imax];

	double ** dx2;
	dx2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dx2[j] = new double[imax];

	double ** ax2;
	ax2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ax2[j] = new double[imax];


	double ** my2;
	my2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my2[j] = new double[imax];

	double ** dy2;
	dy2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dy2[j] = new double[imax];

	double ** ay2;
	ay2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ay2[j] = new double[imax];

	double ** my;
	my = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my[j] = new double[imax];

	double ** mx;
	mx = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mx[j] = new double[imax];


	double ** dy3;
	dy3 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dy3[j] = new double[imax];

	double ** ay3;
	ay3 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ay3[j] = new double[imax];


	double ** dx3;
	dx3 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dx3[j] = new double[imax];

	double ** ax3;
	ax3 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		ax3[j] = new double[imax];

	double ** T;
	T = new double*[jmax];
	for (j = 0; j < jmax; j++)
		T[j] = new double[imax];


	double ** Div;
	Div = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Div[j] = new double[imax];










	double ** A;
	A = new double*[jmax];
	for (j = 0; j < jmax; j++)
		A[j] = new double[imax];

	double ** D;
	D = new double*[jmax];
	for (j = 0; j < jmax; j++)
		D[j] = new double[imax];

	double ** S;
	S = new double*[jmax];
	for (j = 0; j < jmax; j++)
		S[j] = new double[imax];

	double ** p;
	p = new double*[jmax];
	for (j = 0; j < jmax; j++)
		p[j] = new double[imax];


	double ** u_star;
	u_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u_star[j] = new double[imax];

	
		double ** u_old;
		u_old = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u_old[j] = new double[imax];

	double ** v_old;
	v_old = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v_old[j] = new double[imax];


	double ** v_star;
	v_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v_star[j] = new double[imax];

	double ** p_prime;
	p_prime = new double*[jmax];
	for (j = 0; j < jmax; j++)
		p_prime[j] = new double[imax];

	double ** T_old;
	T_old = new double*[jmax];
	for (j = 0; j < jmax; j++)
		T_old[j] = new double[imax];

	double ** v;
	v = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v[j] = new double[imax];

	double ** u;
	u = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u[j] = new double[imax];







	// Coordinates of face centers
	double *x = new double[imax];
	x[0] = 0;
	double C1 = (L1 / (imax - 2));
	for (i = 1; i < imax - 1; i++)
		x[i] = C1*i;
	double *y = new double[jmax];
	y[0] = 0;
	C1 = (L2 / (jmax - 2));
	for (j = 1; j < jmax - 1; j++)
		y[j] = C1*j;



	// Coordinates of Cell Centers; NEEDED FOR PLOTTING
	double *xc = new double[imax];
	xc[0] = 0;
	xc[imax - 1] = L1;
	for (i = 1; i < imax - 1; i++)
		xc[i] = (x[i] + x[i - 1]) / 2.0;
	double *yc = new double[imax];
	yc[0] = 0;
	yc[imax - 1] = L2;
	for (j = 1; j < jmax - 1; j++)
		yc[j] = (y[j] + y[j - 1]) / 2.0;



	//time step
	double Dt = 0.9;
//	double Dt_advec = 1 / ((abs(u) / Dx) + (abs(v) / Dy));
//	double Dt_cond = (1 / 2 * alpha) / (1 / (Dx*Dx) + 1 / (Dy*Dy));
//	double Dt = min(Dt_advec, Dt_cond);



	// velocity and pressure correction BC


	for (j = 0; j < jmax; j++)
	{
		u[j][0] = u_vel_west;
		u[j][imax - 1] = u_vel_east;
		v[j][0] = v_vel_west;
		u[j][imax - 1] = u_vel_east;
		p_prime[j][0] = p_prime[j][1];
		p_prime[j][imax] = p_prime[j][imax - 1];


		
	}

	

	for (i = 0; i < imax; i++)
	{
		u[0][i] = u_vel_south;
		u[jmax - 1][i] = u_vel_north;
		v[0][i] = v_vel_south;
		v[jmax - 1][i] = v_vel_north;
		p_prime[jmax][i] = p_prime[jmax - 1][i];
		p_prime[0][i] = p_prime[1][i];



	}

	





	
	//calc of x-mom

	for (j = 1; j < jmax - 2; j++)
	{
		for (i = 0; i < imax - 1; i++)
		{
			mx1[j][i] = rho*(u[j][i] + u[j][i + 1]) / 2;
			ax1[j][i] = max(mx1[j][i], 0.0)*u[j][i] + min(mx1[j][i], 0.0)*u[j][i + 1];
			dx1[j][i] = myu*(u[j][i + 1] - u[j][i]) / (xc[i + 1] - xc[i]);
		}

	}

	for (j = 0; j < jmax - 1; j++)
	{
		for (i = 1; i < imax - 2;i++)
		{
			my1[j][i] = rho*(v[j][i] + v[j][i + 1]) / 2;
			ay1[j][i] = max(my1[j][i], 0.0)*u[j][i] + min(my1[j][i], 0.0)*u[j][i + 1];
			dy1[j][i] = myu*(u[j][i + 1] - u[j][i]) / (yc[i + 1] - yc[i]);

		}
	}

	for (j = 1; j < jmax - 1; j++)
	{
		for (i = 1; i < imax - 2; i++)
		{

			A[j][i] = (ax1[j][i] - ax1[j][i - 1])*Dy + (ay1[j][i] - ay1[j - 1][i])*Dx;
			D[j][i] = (dx1[j][i] - dx1[j][i - 1])*Dy + (dy1[j][i] - dy1[j - 1][i])*Dx;
			S[j][i] = (p[j][i] - p[j][i + 1])*Dy;
			u_star[j][i] = u_old[j][i] + (Dt / (rho*c*DV))*(D[j][i] - A[j][i] + S[j][i]);
		}

	}


	//calc of y-mom
	for (j = 1; j < jmax - 2; j++)
	{
		for (i = 0; i < imax - 1; i++)
		{
			mx2[j][i] = rho*(u[j][i] + u[j + 1][i]) / 2;
			ax2[j][i] = max(mx2[j][i], 0.0)*v[j][i] - min(mx2[j][i], 0.0)*v[j][i + 1];
			dx2[j][i] = myu*(v[j][i + 1] - v[j][i]) / (xc[i + 1] - xc[i]);


		}
	}


	for (j = 0; j < jmax - 2; j++)
	{
		for (i = 1; i < imax - 1; i++)
		{
			my2[j][i] = rho*(v[j][i] + v[j + 1][i]) / 2;
			ay2[j][i] = max(my2[j][i], 0.0)*v[j][i] - min(my2[j][i], 0.0)*v[j][i + 1];
			dy2[j][i] = myu*(v[j + 1][i] - v[j][i]) / (yc[j + 1] - yc[j]);
			
		}
	}

	for (j = 1; j < jmax - 2; j++)
	{
		for (i = 1; i < imax - 1; i++)
		{


			A[j][i] = (ax2[j][i] - ax2[j][i - 1])*Dy + (ay2[j][i] - ay2[j - 1][i])*Dx;
			D[j][i] = (dx2[j][i] - dx2[j][i - 1])*Dy + (dy2[j][i] - dy2[j - 1][i])*Dx;
			S[j][i] = (p[j][i] - p[j+1][i])*Dy;
			v_star[j][i] = v_old[j][i] + (Dt / (rho*c*DV))*(D[j][i] - A[j][i] + S[j][i]);
		}
	}


	//corrector step

	for (j = 1; j < jmax - 1; j++)
	{

		for (i = 1; i < imax - 1; i++)
		{

			Div[j][i] = (u_star[j][i] - u_star[j][i - 1])*Dy + (v_star[j][i] - v_star[j][i])*Dx;
		}
	}

	for (j = 1; j < jmax; j++)
	{
		for (i = 1; i < imax; i++)
		{


			if (Div[j][i] > conver_crit)  //check converg crit
			{

				//calc of temp

				for (j = 1; j < jmax - 1; j++)
				{
					for (i = 1; i < imax - 1; i++)
					{
						mx[j][i] = rho*u[j][i];
						ax3[j][i] = max(mx[j][i], 0.0)*T[j][i] - min(mx[j][i], 0.0)*T[j][i + 1];
						dx3[j][i] = k*(T[j][i + 1] - T[j][i]) / (xc[j + 1] - xc[j]);


					}
				}

				for (j = 0; j < jmax - 1; j++)
				{
					for (i = 1; i < imax - 1; i++)
					{
						my[j][i] = rho*v[j][i];
						ay3[j][i] = max(my[j][i], 0.0)*T[j][i] - min(my[j][i], 0.0)*T[j + 1][i];
						dy3[j][i] = k*(T[j + 1][i] - T[j][i]) / (yc[i + 1] - yc[i]);



					}
				}

				for (j = 1; j < jmax - 1; j++)
				{
					for (i = 1; i < imax - 1; i++)
					{
						A[j][i] = (ax3[j][i] - ax3[j][i - 1])*Dy + (ay3[j][i] - ay3[j - 1][i])*Dx;
						D[j][i] = (dx3[j][i] - dx3[j][i - 1])*Dy + (dy3[j][i] - dy3[j - 1][i])*Dx;
						S[j][i] = DV*q_gen;
						T[j][i] = T_old[j][i] + (Dt / (c*DV*rho))*(D[j][i] - A[j][i] + S[j][i]);

					}
				}
			}

			else
			{

				continue;
			}

			for (j = 1; j < jmax - 1; j++)
			{
				for (i = 1; i < imax - 1; i++)
				{
					p_prime[j][i] = (aE*p_prime[j][i + 1] + aW*p_prime[j][i - 1] + aN*p_prime[j + 1][i] + aS*p_prime[j + 1][i] - Div[j][i]) / aP;
					p[j][i] = p[j][i] + p_prime[j][i];
				}
			}


			for (j = 1; j < jmax - 1; j++)
			{

				for (i = 1; i < imax - 2; i++)
				{
					u_star[j][i] = u_star[j][i] + (Dt / (rho*DV))*(p_prime[j][i] - p_prime[j][i + 1])*Dy;
				}
			}

			for (j = 1; j < jmax - 2; j++)
			{
				for (i = 1; i < imax - 1; i++)
				{

					v_star[j][i] = v_star[j][i] + (Dt / (rho*DV))*(p_prime[j][i] - p_prime[j + 1][i])*Dx;

				}

			}

		}
	}

}






int main()

{

	void navierstokes_semi_explicit();
	



}


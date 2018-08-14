#include "Rheaders.h"
#include "directions.h"

using namespace std;

void runidir(double *u, double &theta, double &phi)
{
    phi=2.0*PI*runif(0.0,1.0);
    theta=acos(2.0*runif(0.0,1.0)-1.0);

    u[0]=cos(phi)*sin(theta);
    u[1]=sin(phi)*sin(theta);
    u[2]=cos(theta);
}

void vectorRotation(double *u, double *w, double *mu)
{
	double sphi, cphi, sth;
	if(mu[2]*mu[2]<1)
	{
		sth = sqrt(1.0-mu[2]*mu[2]);
		sphi = mu[1]/sth;
		cphi = mu[0]/sth;
		u[0] = mu[2]*cphi*w[0]-sphi*w[1]+sth*cphi*w[2];
		u[1] = mu[2]*sphi*w[0]+cphi*w[1]+sth*sphi*w[2];
		u[2] = -sth*w[0]+mu[2]*w[2];
	}
	else
	{
		if(mu[2]>0) {u[0]=w[0]; u[1]=w[1]; u[2]=w[2];}
		if(mu[2]<0) {u[0]=-w[0]; u[1]=-w[1]; u[2]=-w[2];}
	}

}

void rVonMisesFisher(double *u, double *mu, double kappa, double &theta, double &phi)
{
    double v[2], s, w[3];

    phi = 2.0*PI*runif(0.0,1.0);
    v[0] = cos(phi);
    v[1] = sin(phi);
    w[2] = log(exp(-kappa)+(exp(kappa)-exp(-kappa))*runif(0.0,1.0))/kappa;
    s = sqrt(1-w[2]*w[2]);
    w[0] = s*v[0];
    w[1] = s*v[1];

    vectorRotation(u,w,mu);
    // TODO: get theta from u
    theta = 0.0;
}

void rOhserSchladitz(double *u, double *mu, double kappa, double &theta, double &phi)
{
    double q,w[3];
    phi = 2.0*M_PI*runif(0.0,1.0);  // 0.0; /* DEBUG */
    q = runif(0.0,1.0);
    theta = acos((1-2*q)/sqrt(kappa*kappa-(1-2*q)*(1-2*q)*(kappa*kappa-1)));

    w[0] = sin(theta)*cos(phi);
    w[1] = sin(theta)*sin(phi);
    w[2] = cos(theta);

    vectorRotation(u,w,mu);
}

double scalarProduct(double *x, double *y)
{
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

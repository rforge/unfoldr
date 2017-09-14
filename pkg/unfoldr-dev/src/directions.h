/**
 * @file directions.h
 * @date 2013-11-10
 *
 * @brief random orientation directions
 *
 * @author: F. Ballani
 */

#ifndef DIRECTIONS_H_
#define DIRECTIONS_H_

void runidir(double *u, double &theta, double &phi);
void rVonMisesFisher(double *u, double *mu, double kappa, double &phi);
void rOhserSchladitz(double *u, double *mu, double kappa, double &theta, double &phi);
void vectorRotation(double *u, double *w, double *mu);
double scalarProduct(double *x, double *y);

#endif /* DIRECTIONS_H_ */

#include <math.h>
#include "simpson.hpp"

double DoubleSimpson(double a, double b, double (*f) (double x), double* v)
{
 	double SAB = 0, SAC_SCB = 0, mid, erro, fLo, fHi, fMid;

	mid = (a + b)/2.0;

	fLo = f(a);
	fMid = f(mid);
	fHi = f(b);

	SAB += ((b-a)/6.0) * (fLo + 4*fMid + fHi);

	SAC_SCB += ((mid-a)/6.0) * (fLo + 4*f((a+mid)/2.0) + fMid);
	SAC_SCB += ((b-mid)/6.0) * (fMid + 4*f((mid+b)/2.0) + fHi);

	erro = fabs(SAB - SAC_SCB)/15.0;

	*v = SAC_SCB + (erro);

	return erro;
}

double AdaptiveSimpson(double a, double b, double (*f) (double x), double tol)
{
	double tmp, val;
	tmp = DoubleSimpson(a, b, f, &val);
	if(tmp < tol) return val;
	return AdaptiveSimpson(a, (b+a)/2.0, f, tol/2.0) + AdaptiveSimpson((b+a)/2.0, b, f, tol/2.0);
}

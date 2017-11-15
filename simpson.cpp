#include <cmath>
#include <algorithm>
#include <vector>
#include "simpson.hpp"

using namespace std;

constexpr int Nx = 256;
constexpr int Ny = 256;
constexpr int Nz = 99;

static bool checkInt(double x)
{
	if(fabs(x - floor(x)) < 1e-6) return true;
	return false;
}

static double interpolate(int i, double j, int k, vector<double>& vet)
{
	int lo = (int)floor(j);
	int hi = lo + 1;
	return (vet[k*Nx*Ny + lo*Nx + i] + vet[k*Nx*Ny + hi*Nx + i])/2.0;
}

double simpsonMarretado(int i, int k, double a, double b, vector<double>& vet, double h)
{
	// vamos fazer simpson em intervalos de tamanho = h
	double ans = 0, start, end, mid, fStart, fMid, fEnd;
	int lo, hi;
	int steps = ceil((b - a)/h);
	for(int step = 0; step < steps; ++step)
	{
		start = a + (step)*h;
		end = min(b, start + h);
		mid = (start + end)/2.0;

		lo = (int) floor(start);
		if(!checkInt(start))
		{
			hi = lo+1;
			fStart = (vet[k*Nx*Ny + lo*Nx + i] + vet[k*Nx*Ny + hi*Nx + i])/2.0;
		}
		else
		{
			fStart = vet[k*Nx*Ny + lo*Nx + i];
		}

		lo = (int) floor(mid);
		if(!checkInt(mid))
		{
			hi = lo+1;
			fMid = (vet[k*Nx*Ny + lo*Nx + i] + vet[k*Nx*Ny + hi*Nx + i])/2.0;
		}
		else
		{
			fMid = vet[k*Nx*Ny + lo*Nx + i];
		}

		lo = (int) floor(end);
		if(!checkInt(end))
		{
			hi = lo+1;
			fEnd = (vet[k*Nx*Ny + lo*Nx + i] + vet[k*Nx*Ny + hi*Nx + i])/2.0;
		}
		else
		{
			fEnd = vet[k*Nx*Ny + lo*Nx + i];
		}
		resp += (h/6.0) * (fStart + 4*fMid + fEnd);
	}

	return resp;
}



double DoubleSimpson(int a, int b, vector<double>& vet, double* v)
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

#include <bits/stdc++.h>

using namespace std;

struct Simpson{
	vector<double> image;

	Simpson(vector<double>& vet);

	double tau(double s);

	double interpolate(int i, int k, double mid);

	double simpsonFilho(int i, int k, double s, int h);

	double simpsonPai(int i, int k, int b, int h);

}
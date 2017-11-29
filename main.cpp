#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

constexpr int Nx = 256;
constexpr int Ny = 256;

constexpr double TOL = 1e-8;
constexpr double minTol = 1e-7;

int totalDeChamadas = 0;

inline bool checkInt(double x)
{
	if(fabs(x - floor(x)) < 1e-4) return true;
	return false;
}

inline int cmpDouble(double& a, double& b){
	return (a + TOL > b) ? ((b + TOL > a) ? 0 : 1) : -1;
}

inline int idx(int a, int b, int c){
	return a*Nx*Nx + b*Nx + c;
}

inline void updateTol(double& x){
	x = (x < minTol)? minTol : x;
}

class Simpson{

public:
	vector<double> img;
	
	Simpson(vector<double>& vet){
		img = vet;
	}

	inline double tau(double s)
	{
		return (s < 0.30) ? 0.0 : (0.05 * (s - 0.3));
	}

	double interpolate(int i, int k, double mid)
	{
		int lo = (int) floor(mid);
		int hi = lo+1;
		double w1 = 1.0 - (mid - lo);
		double w2 = 1.0 - w1;
		return (w1 * tau(img[idx(k, lo, i)]) + w2 * tau(img[idx(k, hi, i)]) );
	}

	double simpsonFilhoAdaptativo(int i, int k, double a, double b, double tol)
	{
		updateTol(tol);

		totalDeChamadas++;

		double c = (a + b)/2.0;
		double fA, fB, fC, fAC, fBC;

		double h = b - a;

		double SAB, SAC_SCB;

		if(!checkInt(a))
		{
			//printf("A\n");
			fA = interpolate(i, k, a);
		}
		else
		{
			//printf("B\n");
			fA = tau(img[idx(k, (int)a, i)]);
		}

		if(!checkInt(b))
		{
			//printf("A1\n");
			fB = interpolate(i, k, b);
		}
		else
		{
			//printf("B1\n");
			fB = tau(img[idx(k, (int)b, i)]);
		}

		if(!checkInt(c))
		{
			//printf("A2\n");
			fC = interpolate(i, k, c);
		}
		else
		{
			//printf("B2\n");
			fC = tau(img[idx(k, (int)c, i)]);
		}

		SAB = (h/6.0) * (fA + 4*fC + fB);

		double midAC, midCB;

		midAC = (a + c)/2.0;
		midCB = (c + b)/2.0;

		if(!checkInt(midAC))
		{
			fAC = interpolate(i, k, midAC);
		}
		else
		{
			fAC = tau(img[idx(k, (int)midAC, i)]);
		}

		if(!checkInt(midCB))
		{
			fBC = interpolate(i, k, midCB);
		}
		else
		{
			fBC = tau(img[idx(k, (int)midCB, i)]);
		}

		//cout << "XEGUEI AQUI MIZERAVI" << endl;
		SAC_SCB = (h/12.0) * (fA + 4*fAC + fC);
		SAC_SCB += (h/12.0) * (fC + 4*fBC + fB);

		double error = fabs(SAB - SAC_SCB)/15.0;

		if(cmpDouble(error, tol) <= 0)
		{
			return SAC_SCB;
		}
		else
		{
			return simpsonFilhoAdaptativo(i, k, a, c, (tol/2.0)) + simpsonFilhoAdaptativo(i, k, c, b, (tol/2.0));
		}

	}

	double adaptativeSimpson(int i, int k, double a, double b, double tol)
	{
		updateTol(tol);

		totalDeChamadas++;

		double mid = (a+b)/2.0;
		double lowMid = (a+mid)/2.0;
		double highMid = (mid + b)/2.0;
		double h = b - a;

		double resultSimpsonFilhoStart = -simpsonFilhoAdaptativo(i, k, 0.0, a, tol);
		double resultSimpsonFilhoMid = -simpsonFilhoAdaptativo(i, k, 0.0, mid, tol);
		double resultSimpsonFilhoEnd = -simpsonFilhoAdaptativo(i, k, 0.0, b, tol);
		double resultSimpsonFilhoLowMid = -simpsonFilhoAdaptativo(i, k, 0.0, lowMid, tol);
		double resultSimpsonFilhoHighMid = -simpsonFilhoAdaptativo(i, k, 0.0, highMid, tol);

		double fA, fB, fMid, fLowMid, fHighMid;


		if(!checkInt(a))
		{
			fA = interpolate(i, k, a);
		}
		else
		{
			fA = tau(img[idx(k, (int)a, i)]);
		}
		fA *= exp(resultSimpsonFilhoStart);


		if(!checkInt(b))
		{
			fB = interpolate(i, k, b);
		}
		else
		{
			fB = tau(img[idx(k, (int)b, i)]);
		}
		fB *= exp(resultSimpsonFilhoEnd);


		if(!checkInt(mid))
		{
			fMid = interpolate(i, k, mid);
		}
		else
		{
			fMid = tau(img[idx(k,(int)mid, i)]);
		}
		fMid *= exp(resultSimpsonFilhoMid);


		double SAB = (h/6.0) * (fA + 4*fMid + fB);

		if(!checkInt(lowMid))
		{
			fLowMid = interpolate(i, k, lowMid);
		}
		else
		{
			fLowMid = tau(img[idx(k, (int)lowMid, i)]);
		}
		fLowMid *= exp(resultSimpsonFilhoLowMid);


		if(!checkInt(highMid))
		{
			fHighMid = interpolate(i, k, highMid);
		}
		else
		{
			fHighMid = tau(img[idx(k, (int)highMid, i)]);
		}
		fHighMid *= exp(resultSimpsonFilhoHighMid);


		double SAC_SCB = (h/12.0) * (fA + 4*fLowMid + fMid);
		SAC_SCB += (h/12.0) * (fMid + 4*fHighMid + fB);

		double error = fabs(SAC_SCB - SAB)/15.0;

		if(cmpDouble(error, tol) <= 0)
		{
			return SAC_SCB;
		}
		return adaptativeSimpson(i, k, a, mid, tol/2.0) + adaptativeSimpson(i, k, mid, b, tol/2.0);
	}

	// Estamos colocando h como step inteiro inicialmente
	double simpsonFilho(int i, int k, double s, double h)
	{
		totalDeChamadas++;
		//cout << "entrei no simpson filho" << endl;
		double ans = 0;
		double fStart, fMid, fEnd, mid, start, lo, hi, end;
		//int lo, hi, start, end;
		int steps = ceil(s/h);
		//cout << "s = " << s << endl;
		for(int step = 0; step < steps; ++step)
		{

			//cout << "estou no step = " << step << endl;

			start = step*h;
			end = min(s, start + h);
			mid = (start + end) / 2.0;

			//if(step == 0) fStart = tau(img[idx(k, start, i)]);
			//else fStart = fEnd;
			
			if(step != 0) fStart = fEnd;
			else if(!checkInt(start)){
				fStart = interpolate(i,k,start);
			}
			else fStart = tau(img[idx(k, (int)start, i)]);
			//cout << "cheguei aqui 1 " << endl;		

			//cout << "cheguei aqui 2 " << endl;
			if(!checkInt(mid))
			{
				fMid = interpolate(i, k, mid);
			}
			else
			{
				fMid = tau(img[idx(k, (int)mid, i)]);
			}

			if(!checkInt(end))
			{
				fEnd = interpolate(i, k, end); 
			}
			else fEnd = tau(img[idx(k, (int)end, i)]);
			

			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
			//cout << "sai do SimpsonFilho" << endl;
		}
		return ans;
	}

	double simpsonPai(int i, int k, double b, double h)
	{
		totalDeChamadas++;
		double ans = 0;
		double fStart, fMid, fEnd, mid;
		double lo, hi, start, end;
		
		int steps = ceil(b/(double)h);

		//cout << "numero de steps = " << steps << endl;
		for(int step = 0; step < steps; ++step)
		{
			start = step*h;
			end = min(b, start + h);
			mid = (start + end) / 2.0;

			double resultSimpsonFilhoStart = -simpsonFilho(i, k, start, h);
			double resultSimpsonFilhoMid = -simpsonFilho(i, k, mid, h);
			double resultSimpsonFilhoEnd = -simpsonFilho(i, k, end, h);

			if(step == 0) fStart = tau(img[idx(k, start, i)]) * exp(resultSimpsonFilhoStart);
			else fStart = fEnd;

			if(step != 0)
			{
				fStart = fEnd;
			}
			else if(!checkInt(start)){
				fStart = interpolate(i, k, start) * exp(resultSimpsonFilhoStart);
			}
			else
			{
				fStart = tau(img[idx(k, (int)start, i)]) *  exp(resultSimpsonFilhoStart);
			}


			if(!checkInt(mid))
			{
				fMid = interpolate(i, k, mid) * exp(resultSimpsonFilhoMid);
			}
			else
			{
				fMid = tau(img[idx(k, (int)mid, i)]) * exp(resultSimpsonFilhoMid);
			}


			if(!checkInt(end))
			{
				fEnd = interpolate(i,k,end) * exp(resultSimpsonFilhoEnd);
			}
			else
			{
				fEnd = tau(img[idx(k, (int) end, i)]) * exp(resultSimpsonFilhoEnd);
			}
			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
			//cout << "terminei o step" << endl;
		}
		return ans;
	}

};

vector<unsigned char> readData(const string& fileName)
{	
	ifstream file(fileName, ios::in | ios::binary);
	vector<unsigned char> contents((istreambuf_iterator<char>(file)),
		istreambuf_iterator<char>());
	return contents;
}

int main()
{
	vector<unsigned char> data = readData("head-8bit.raw");
	vector<double> cpy;
	ofstream out;
	out.open("saidaAdaptativo.pgm", ofstream::out);

	for(const unsigned int& x : data) cpy.push_back((double)x/255.0); // vetor normalizado [0,1]
	
	Simpson teste(cpy);

	unsigned char mx = 0;
	out << "P2" << endl;
	out << "128 99" << endl;

	vector<vector<unsigned char> > imagem(128, vector<unsigned char>(99, 0.0));
	//vector<vector<unsigned char> > imagem2(128, vector<unsigned char>(128, 0.0));


	for(int i = 0; i < 128; ++i)
	{
		for(int k = 0; k < 99; ++k)
		{
			// double t1 = teste.simpsonPai(2*i, k, 255, 4.5);
			// double t2 = teste.simpsonPai(2*i + 1, k, 255, 4.5); 
			// imagem[i][k] = (unsigned char) round(((t1 + t2)/2.0) * 255.0);
			// mx = max(imagem[i][k], mx);

			// Descomentar para chamar Simpson Adaptativo !
			double t1 = teste.adaptativeSimpson(2*i, k, 0, 255, 1e-5);
			double t2 = teste.adaptativeSimpson(2*i + 1, k, 0, 255, 1e-5);
			imagem[i][k] = (unsigned char) round(((t1 + t2)/2.0) * 255.0);
			mx = max(imagem[i][k], mx);
		}
	}

	out << (unsigned int) mx << endl;

	for(int coluna = 0; coluna < 99; ++coluna)
	{
		for(int linha = 0; linha < 128; ++linha)
		{
			out << (unsigned int) imagem[linha][coluna] << " ";
		}
		out << endl;
	}

	cout << "Fiz um total de " << totalDeChamadas << " chamadas de simpson" << endl;
	
	out.close();

	return 0;
}

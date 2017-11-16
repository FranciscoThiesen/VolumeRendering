#include <bits/stdc++.h>
#include <cmath>

using namespace std;

constexpr int Nx = 256;
constexpr int Ny = 256;

bool checkInt(double x)
{
	if(fabs(x - floor(x)) < 1e-6) return true;
	return false;
}

class Simpson{

public:
	vector<double> img;
	
	Simpson(vector<double>& vet){
		img = vet;
	}

	double tau(double s)
	{
		return (s < 0.30) ? 0.0 : (0.05 * (s - 0.3));
	}


	double interpolate(int i, int k, double mid)
	{
		int lo = (int) floor(mid);
		int hi = lo+1;
		return (tau(img[k*Nx*Nx + lo*Nx + i]) + tau(img[k*Nx*Nx + hi*Nx + i]))/2.0;
	}

	// Estamos colocando h como step inteiro inicialmente
	double simpsonFilho(int i, int k, double s, int h)
	{
		//cout << "entrei no simpson filho" << endl;
		double ans = 0;
		double fStart, fMid, fEnd, mid;
		int lo, hi, start, end;
		int steps = ceil(s/h);
		//cout << "s = " << s << endl;
		for(int step = 0; step < steps; ++step)
		{

			//cout << "estou no step = " << step << endl;

			start = step*h;
			end = min((int)s, start + h);
			mid = (start + end) / 2.0;

			if(step == 0) fStart = tau(img[k*Nx*Ny + start*Nx + i]);
			else fStart = fEnd;
			//cout << "cheguei aqui 1 " << endl;		

			//cout << "cheguei aqui 2 " << endl;

			if(!checkInt(mid))
			{
				fMid = interpolate(i, k, mid);
			}
			else
			{
				int aux = (int)mid;
				fMid = tau(img[k*Nx*Nx + aux*Nx + i]);
			}

			//cout << "cheguei aqui 3 " << endl;

			fEnd = tau(img[k*Nx*Ny + end*Nx + i]);

			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
			//cout << "sai do SimpsonFilho" << endl;
		}
		return ans;
	}

	double simpsonPai(int i, int k, int b, int h)
	{
		double ans = 0;
		double fStart, fMid, fEnd, mid;
		int lo, hi, start, end, aux;
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

			if(step == 0) fStart = tau(img[k*Nx*Ny + start*Nx + i]) * exp(resultSimpsonFilhoStart);
			else fStart = fEnd;

			if(!checkInt(mid))
			{
				double prev = interpolate(i, k, mid);
				fMid = prev * exp(resultSimpsonFilhoMid);
			}
			else
			{
				aux = (int)mid;
				fMid = tau(img[k*Nx*Nx + aux*Nx + i]) * exp(resultSimpsonFilhoMid);
			}

			fEnd = tau(img[k*Nx*Ny + end*Nx + i]) * exp(resultSimpsonFilhoEnd);

			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
			//cout << "terminei o step" << endl;
		}
		return ans;
	}
};

// Function for reading CT scans as binary mode file
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
	out.open("saida.pgm", ofstream::out);

	for(const unsigned int& x : data) cpy.push_back((double)x/255.0); // vetor normalizado [0,1]
	
	Simpson teste(cpy);

	unsigned char mx = 0;
	out << "P2" << endl;
	out << "128 99" << endl;

	vector<vector<unsigned char> > imagem(128, vector<unsigned char>(128, 0.0));

	for(int i = 0; i < 128; ++i){
		for(int k = 0; k < 99; ++k){
			double t1 = teste.simpsonPai(2*i, k, 255, 4);
			double t2 = teste.simpsonPai(2*i + 1, k, 255, 4);
			imagem[i][k] = (unsigned char) round(((t1 + t2)/2.0) * 255.0);
			mx = max(imagem[i][k], mx);
		}
	}

	out << (unsigned int)mx << endl;

	// imagem estava com as dimensoes invertidas, vamos escrever a transposta
	for(int coluna = 0; coluna < 99; ++coluna){
		for(int linha = 0; linha < 128; ++linha){
			out << (unsigned int) imagem[linha][coluna] << " ";
		}
		out << endl;
	}
	out.close();
	return 0;
}

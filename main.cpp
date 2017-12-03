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

//Método que verifica se um valor double pode ser considerado inteiro dada a tolerância 1e-4
//Params: x -> valor double a ser verificado
inline bool isInt(double x)
{
	return (fabs(x - floor(x)) < 1e-4);
}

//Método que compara dois valores double considerando uma tolerância de TOL
//Params: a, b -> valores double a serem comparados
//Return: se a < b -> -1, se a == b -> 0, se a > b -> 1
inline int cmpDouble(double a, double b){
	return (a + TOL > b) ? ((b + TOL > a) ? 0 : 1) : -1;
}

//Método que retorna a posição do vetor onde está o valor associado à uma coordenada da matriz tridimensional
//Params: a, b, c -> coordenadas do ponto na matriz tridimensional na configuração (x,y,z)
inline int idx(int a, int b, int c){
	return a*Nx*Nx + b*Nx + c;
}

//Método que atualiza a tolerância considerando uma tolerância mínima
//Params: x -> referência para o valor atual da tolerância
//Return: se x < minTol -> x = minTol, caso contrário x continua inalterado
inline void updateTol(double& x){
	x = (x < minTol)? minTol : x;
}

class Simpson{

public:
	vector<double> img;
	
	// Construtor da classe Simpsons
	// Params: vet -> vetor normalizado entre 0 e 1
	Simpson(vector<double>& vet){
		img = vet;
	}

	// Método de transferência
	// Params: s -> valor float de 0 a 1 a ser avaliado
	inline double tau(double s)
	{
		return (s < 0.30) ? 0.0 : (0.05 * (s - 0.3));
	}

	// Método que calcula o valor da função tau associada a um ponto que não é inteiro usando a média ponderada
	// dos valores da função tau dos seus vizinhos inteiros.
	// Params: i, k -> posições da matriz tridimensional onde estamos executando a integral
	//		   mid -> valor não inteiro que será usado para determinar os pesos associados a cada vizinho
	double interpolate(int i, int k, double mid)
	{
		int lo = (int) floor(mid);
		int hi = lo+1;
		double w1 = 1.0 - (mid - lo);
		double w2 = 1.0 - w1;
		return (w1 * tau(img[idx(k, lo, i)]) + w2 * tau(img[idx(k, hi, i)]) );
	}

	// Método que executa a integração numérica através do método de Simpson adaptativo da função presente
	// no expoente da expressão para cálculo de intensidade.
	// Params: i, k -> posições da matriz tridimensional onde estamos executando a integral
	//		   a -> limite inferior da chamada corrente
	//		   b -> limite superior da chamada corrente
	//		   tol -> tolerância máxima permitida para aceitação do valor de retorno
	double simpsonAdaptativoAuxiliar(int i, int k, double a, double b, double tol)
	{
		updateTol(tol);

		totalDeChamadas++;

		double c = (a + b)/2.0;
		double fA, fB, fC, fAC, fBC;

		double h = b - a;

		double SAB, SAC_SCB;

		if(!isInt(a))
		{
			fA = interpolate(i, k, a);
		}
		else
		{
			fA = tau(img[idx(k, (int)a, i)]);
		}

		if(!isInt(b))
		{
			fB = interpolate(i, k, b);
		}
		else
		{
			fB = tau(img[idx(k, (int)b, i)]);
		}

		if(!isInt(c))
		{
			fC = interpolate(i, k, c);
		}
		else
		{
			fC = tau(img[idx(k, (int)c, i)]);
		}

		SAB = (h/6.0) * (fA + 4*fC + fB);

		double midAC, midCB;

		midAC = (a + c)/2.0;
		midCB = (c + b)/2.0;

		if(!isInt(midAC))
		{
			fAC = interpolate(i, k, midAC);
		}
		else
		{
			fAC = tau(img[idx(k, (int)midAC, i)]);
		}

		if(!isInt(midCB))
		{
			fBC = interpolate(i, k, midCB);
		}
		else
		{
			fBC = tau(img[idx(k, (int)midCB, i)]);
		}

		SAC_SCB = (h/12.0) * (fA + 4*fAC + fC);
		SAC_SCB += (h/12.0) * (fC + 4*fBC + fB);

		double error = fabs(SAB - SAC_SCB)/15.0;

		if(cmpDouble(error, tol) <= 0)
		{
			return SAC_SCB;
		}
		else
		{
			return simpsonAdaptativoAuxiliar(i, k, a, c, (tol/2.0)) + simpsonAdaptativoAuxiliar(i, k, c, b, (tol/2.0));
		}

	}

	//Método para calcular a integral numéricamente usando o método de Simpson adaptativo.
	//Params: i, k -> posições da matriz tridimensional onde estamos executando a integral
	//		   a -> limite inferior da chamada corrente
	//		   b -> limite superior da chamada corrente
	//		   tol -> tolerância máxima permitida para aceitação do valor de retorno
	double simpsonAdaptativoPrincipal(int i, int k, double a, double b, double tol)
	{
		updateTol(tol);

		totalDeChamadas++;

		double mid = (a+b)/2.0;
		double lowMid = (a+mid)/2.0;
		double highMid = (mid + b)/2.0;
		double h = b - a;

		double resultSimpsonFilhoStart = -simpsonAdaptativoAuxiliar(i, k, 0.0, a, tol);
		double resultSimpsonFilhoMid = -simpsonAdaptativoAuxiliar(i, k, 0.0, mid, tol);
		double resultSimpsonFilhoEnd = -simpsonAdaptativoAuxiliar(i, k, 0.0, b, tol);
		double resultSimpsonFilhoLowMid = -simpsonAdaptativoAuxiliar(i, k, 0.0, lowMid, tol);
		double resultSimpsonFilhoHighMid = -simpsonAdaptativoAuxiliar(i, k, 0.0, highMid, tol);

		double fA, fB, fMid, fLowMid, fHighMid;


		if(!isInt(a))
		{
			fA = interpolate(i, k, a);
		}
		else
		{
			fA = tau(img[idx(k, (int)a, i)]);
		}
		fA *= exp(resultSimpsonFilhoStart);


		if(!isInt(b))
		{
			fB = interpolate(i, k, b);
		}
		else
		{
			fB = tau(img[idx(k, (int)b, i)]);
		}
		fB *= exp(resultSimpsonFilhoEnd);


		if(!isInt(mid))
		{
			fMid = interpolate(i, k, mid);
		}
		else
		{
			fMid = tau(img[idx(k,(int)mid, i)]);
		}
		fMid *= exp(resultSimpsonFilhoMid);


		double SAB = (h/6.0) * (fA + 4*fMid + fB);

		if(!isInt(lowMid))
		{
			fLowMid = interpolate(i, k, lowMid);
		}
		else
		{
			fLowMid = tau(img[idx(k, (int)lowMid, i)]);
		}
		fLowMid *= exp(resultSimpsonFilhoLowMid);


		if(!isInt(highMid))
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
		return simpsonAdaptativoPrincipal(i, k, a, mid, tol/2.0) + simpsonAdaptativoPrincipal(i, k, mid, b, tol/2.0);
	}

	//Método para cálculo da integral presente no expoente da função de intensidade usando o método de Simpson
	//com passo fixo.
	//Params: i, k -> posições da matriz tridimensional onde estamos executando a integral
	//		   s -> limite superior da integral
	//		   h -> passo tomado a cada iteração
	double simpsonAuxiliar(int i, int k, double s, double h)
	{
		totalDeChamadas++;
		double ans = 0;
		double fStart, fMid, fEnd, mid, start, hi, end;
		int steps = ceil(s/h);
		for(int step = 0; step < steps; ++step)
		{
			start = step*h;
			end = min(s, start + h);
			mid = (start + end) / 2.0;
			
			if(step != 0) fStart = fEnd;
			else if(!isInt(start)){
				fStart = interpolate(i,k,start);
			}
			else fStart = tau(img[idx(k, (int)start, i)]);
			
			if(!isInt(mid))
			{
				fMid = interpolate(i, k, mid);
			}
			else
			{
				fMid = tau(img[idx(k, (int)mid, i)]);
			}

			if(!isInt(end))
			{
				fEnd = interpolate(i, k, end); 
			}
			else fEnd = tau(img[idx(k, (int)end, i)]);
			

			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
		}
		return ans;
	}

	//Método que calcula a integral numericamente utilizando método de Simpson com passo fixo
	//Params: i, k -> posições da matriz tridimensional onde estamos executando a integral
	//		   b -> limite superior da integral a ser calculada
	//		   h -> tamanho do passo a ser tomado em cada iteração
	double simpsonPrincipal(int i, int k, double b, double h)
	{
		totalDeChamadas++;
		double ans = 0;
		double fStart, fMid, fEnd, mid;
		double hi, start, end;
		
		int steps = ceil(b/(double)h);

		for(int step = 0; step < steps; ++step)
		{
			start = step*h;
			end = min(b, start + h);
			mid = (start + end) / 2.0;

			double resultSimpsonFilhoStart = -simpsonAuxiliar(i, k, start, h);
			double resultSimpsonFilhoMid = -simpsonAuxiliar(i, k, mid, h);
			double resultSimpsonFilhoEnd = -simpsonAuxiliar(i, k, end, h);

			if(step == 0) fStart = tau(img[idx(k, start, i)]) * exp(resultSimpsonFilhoStart);
			else fStart = fEnd;

			if(step != 0)
			{
				fStart = fEnd;
			}
			else if(!isInt(start)){
				fStart = interpolate(i, k, start) * exp(resultSimpsonFilhoStart);
			}
			else
			{
				fStart = tau(img[idx(k, (int)start, i)]) *  exp(resultSimpsonFilhoStart);
			}


			if(!isInt(mid))
			{
				fMid = interpolate(i, k, mid) * exp(resultSimpsonFilhoMid);
			}
			else
			{
				fMid = tau(img[idx(k, (int)mid, i)]) * exp(resultSimpsonFilhoMid);
			}


			if(!isInt(end))
			{
				fEnd = interpolate(i,k,end) * exp(resultSimpsonFilhoEnd);
			}
			else
			{
				fEnd = tau(img[idx(k, (int) end, i)]) * exp(resultSimpsonFilhoEnd);
			}
			ans += (h/6.0) * (fStart + 4*fMid + fEnd);
		}
		return ans;
	}

};

//Método que lê arquivo de forma binária e produz vetor de unsigned char
//Params: fileName -> nome do arquivo a ser lido
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

	// Normalização do vetor lido para [0,1]
	for(const unsigned int& x : data) cpy.push_back((double)x/255.0); 
	
	Simpson teste(cpy);

	unsigned char mx = 0;
	out << "P2" << endl;
	out << "128 99" << endl;

	vector<vector<unsigned char> > imagem(128, vector<unsigned char>(99, 0.0));

	for(int i = 0; i < 128; ++i)
	{
		for(int k = 0; k < 99; ++k)
		{
			// double t1 = teste.simpsonPrincipal(2*i, k, 255, 4.5);
			// double t2 = teste.simpsonPrincipal(2*i + 1, k, 255, 4.5); 
			// imagem[i][k] = (unsigned char) round(((t1 + t2)/2.0) * 255.0);
			// mx = max(imagem[i][k], mx);

			// Descomentar para chamar Simpson Adaptativo !
			double t1 = teste.simpsonAdaptativoPrincipal(2*i, k, 0, 255, 1e-3);
			double t2 = teste.simpsonAdaptativoPrincipal(2*i + 1, k, 0, 255, 1e-3);
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

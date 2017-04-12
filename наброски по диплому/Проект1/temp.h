#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

static const double eps = 1e-12;
static const double x_0 = 0;
static const double y_0 = 1;
static const double a1 = 3;
static const double a2 = 2;
static const double b1 = 1;
static const double b2 = 2;
static const double d1 = 3;
static const double d2 = 1;
static const double c1 = 1;
static const double c2 = 2;
static const double A1 = 1;
static const double A2 = 2;
static int num_f=0;

double a = 0, b = 0;
//int num_f = 0;
//lab3 begin
double betta;
static const int mode = 1;
//lab3 end

class vect
{
private:
	double x[2];
public:
	inline vect()
	{
		x[0] = x[1] = 0.0;
	}
	inline vect(double z, double y)
	{
		x[0] = z;
		x[1] = y;
	}
	inline double norm()
	{
		return sqrt(x[0] * x[0] + x[1] * x[1]);
	}
	inline double& operator [] (unsigned i)
	{
		return x[i];
	}
	inline vect operator + (vect y)
	{
		return vect(x[0] + y[0], x[1] + y[1]);
	}
	inline vect operator - (vect y)
	{
		return vect(x[0] - y[0], x[1] - y[1]);
	}
	inline vect operator * (double c)
	{
		return vect(x[0] * c, x[1] * c);
	}
	inline vect operator / (double c)
	{
		return vect(x[0] / c, x[1] / c);
	}
	inline double mult(vect y)
	{
		return x[0] * y[0] + x[1] * y[1];
	}
	friend ostream& operator << (ostream& ostream_, const vect& v)
	{
		ostream_.setf(ios::scientific);
		ostream_.precision(16);
		ostream_ << "[ " << v.x[0] << ", " << v.x[1] << " ]";
		return ostream_;
	}
};

double G1(vect xk)
{
	double result;

	result = xk[0] - 3*xk[1] -3;

	return(result);
}

inline double f(class vect x)
{
//	f_calc++;
	//return((1 - x[0])*(1 - x[0]) + 100 * pow((x[1] - x[0] * x[0]), 2));
	//return( pow((x[0]),2) + pow((x[1]),2));
	return(-A1*exp(-pow((x[0] - a1) / b1, 2) - pow(x[1] - c1, 2) / d1) - A2*exp(-pow((x[0] - a2) / b2, 2) - pow((x[1] - c2) / d2, 2)));
}

//G(q)
//q = 3*x + y - 3
inline double G(class vect x, double rk)
{
	if (mode == 1)
	{
		double check = (0.5*(G1(x) + abs(G1(x))));
		return(0.5*(G1(x) + abs(G1(x))));
	}
	else
		return(-1 / G1(x));
}

//Q(x,r)
inline double Q(class vect x, double rk)
{
	double check = f(x);
	check += rk*G(x, rk);
	return(f(x)+ rk*G(x, rk)); // CHECK IT UP. f + r*G OR f + r*r*G // what is r exactly for?
}

void interval( class vect x, double r, int i)
{
	double eps = 1E-12;
	double h = 1E-07;
	vect xp, x1;
	int n = 2;

	int t = 0;
	for (int j = 0; j < n; j++)
		xp[j] = x1[j] = x[j];
	double fp, ff, f;
	f = Q(x1, r); num_f++;
	//num_f++;
	x1[i] -= h;
	fp = Q(x1, r); // слева
	num_f++;
	//num_f++;
	x1[i] = x[i] + h;
	ff = Q(x1, r); // справа
	num_f++;
	//num_f++;
	x1[i] = x[i];
	if ((fp - f > eps) && (ff - f > eps))
	{
		a = x1[i] - h;
		b = x1[i] + h;
	}
	else
	{
		if (f - ff>eps)
			x1[i] += h, t = 1;
		else h = -h, x1[i] += h;

		fp = Q(xp, r); num_f++;
		ff = Q(x1, r); num_f++;

		int j = 0;
		while (fp - ff>eps)
		{
			j++;
			xp[i] = x1[i];
			fp = ff;
			h *= 2;
			x1[i] += h;
			//
			if (mode == 1)
			{
				ff = Q(x1, r); num_f++;
			}
			else
			{
				if (G1(x1) <= 0)
				{
					ff = Q(x1, r);
					num_f++;
				}
				else
					while (G1(x1) > 0 && h > 1e-20)
						x1[i] -= h, h /= 2, x1[i] += h;
			}
			//
		}
		if (t == 0)
			a = x1[i], b = xp[i] - h / 2;
		else
			b = x1[i], a = xp[i] - h / 2;
		for (int j = 0; j< n; j++) x1[j] = x[j];
	}
}

double dihotomia(int i, double r, vect x)
{
	int n = 2;
	vect xp;
	double eps = 1E-07;
	double d = eps / 4;
	interval(x ,r, i);
	for (int j = 0; j < n; j++)
		xp[j] = x[j];
	double ak = a, bk = b;
	for (int j = 1; fabs(bk - ak) > eps; j++)
	{
		xp[i] = (ak + bk - d) / 2;
		x[i] = (ak + bk + d) / 2;
		if (Q(xp, r) - Q(x, r)<eps) bk = x[i];
		else ak = xp[i];
		num_f += 2;
	}
	x[i] = (ak + bk) / 2;
	return x[i];
}

/**
bool stop_condition(class vect x1, class vect x, double rk)
{
	double check1 = Q(x1, rk);
	double check2 = Q(x, rk);
	if (fabs(Q(x1, rk) - Q(x, rk)) - eps <0.0001)
		return true;
	return false;
}
*/

void main()
{

	if (mode == 1)
		betta = 2;
	else
		betta = 0.5;

	class vect x, xold, s[2], a[2], b;
	double lambda1, lambda2;
	int iter = 0;
//	f_calc = 0;
	double rk = 1;
	x[0] = x_0;
	x[1] = y_0;

	s[0][0] = s[1][1] = 1.0;
	s[0][1] = s[1][0] = 0.0;
	double check;
	check = rk*G(x, rk);
	do
		{
		double eps = 1E-11;
		int k;
		double f = Q(x, rk); num_f += 2;
		double fk = f + 3;
		for (k = 0; k < 1000 && fabs(f - fk)>eps; k++)
		{
			//betta = 2;
			x[0] = dihotomia(0, rk, x);
			x[1] = dihotomia(1, rk, x);
		}
		iter++;
		
		//double check;
		check = rk*G(x, rk);
		rk = rk*betta;

			cout << num_f << '\n' << iter << '\n' << x << '\t' << endl;
	} while (abs(check) > eps);
	cout << f(x);
	system("pause");

}


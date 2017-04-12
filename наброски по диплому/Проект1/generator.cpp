#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Windows.h>
#include <conio.h>

#include <random>

using namespace std;
////////////////////////////////////////////////////
struct point
{
	double x,y;
};
////////////////////////////////////////////////////
// ��������� "����������" ������
vector<double> gen_norm_error(int size)
{
	// random device class instance, source of 'true' randomness for initializing random seed
	std::random_device rd;
	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen(rd());

	vector<double> res;
	res.reserve(size);

	for (int i = 0; i < size; i++) // ��� ���� ��������� �����
	{
		std::normal_distribution<double> d(0, 1); // ����������� �������������
		
		res[i] = d(gen);
	}
	return res;
}
////////////////////////////////////////////////////
// ��������� ������� ������
class data_generation
{
private:
	vector<std::vector<point> > grid; // ����� 
	vector<double> D, q_barriers; // ������ D � ������;
	double x_true, y_true;
	int t;
public:
	double calc_corrupted_d(double x,double y, double error)
	{
		return sqrt( pow(x - x_true, 2) + pow(y - y_true,2) ) + error;
	}
	void Calc_D(vector<double> d_array)//������, ��� ����������� ������� ������� ���������� �������.
	{
		int dim = q_barriers.size();
		int grid_size = grid[0].size * grid[0].size; // ������� ����!!!
		
		for (int i = 0; i < grid_size; i++)
		{
			if (d_array[i] < q_barriers[0]) // ��������� ������� ����������� � ���������� ��������
			{
				D[i] = 0;
				break;
			}
			if (d_array[i] >= q_barriers[dim]) //!!! ����� ��� ��������� ������������ �����
			{
				D[i] = dim;
				break;
			}
			for (int j = 1, int j_next = 2; j < dim - 1; j++, j_next++)
			{
				
			}
		}
	}
	data_generation(double a, double b, int n, int size)//size - ������ �����������. �� [a;b] ������������� n �����. ��� ������ � t(����� ������)?
	{
		//��������� ������
		grid.reserve(n*n);
		q_barriers.reserve(size);
		D.reserve(size + 1);

		//���-�� ������ ��������� �������		
		x_true = (b - a) / 2;
		y_true = (b - a) / 2;
		//
		//������ ����������� - ��������, ����� �� ����(��� �������� �������, ����� ����� ����������� ������������� �������).
		q_barriers[0] = (b - a) / 2;
		//
		//������ ����� xi yi � ����� ����, ����� ��� �������������.
		double step = (b - a) / n;

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				grid[i][j].x = a + i*step;
				grid[i][j].y = a + j*step;
			}
		}
		//
		// ���������� N^2 ������(��� ������� ���� �����)
		vector<double> err = gen_norm_error(n*n);
		//
		// ������ ���� d - ���������� �� ���� �� �������
		vector<double> d;
		d.reserve(n*n);

		for(int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				d[i + j] = calc_corrupted_d(grid[i][j].x,grid[i][j].y,err[i + j]);
			}
		}
		//
		//������ ���� D
		Calc_D(d);
		//
	}
};
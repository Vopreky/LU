#include <iostream>
#include <chrono>
#include <ctime>
#include <cblas.h>//BLAS!
using namespace std;

class Matrix
{
	double * arr;
	int row_count, col_count;

Matrix(const Matrix & M) { throw 0; }
Matrix & operator=(const Matrix & M) { throw 0; }
public:
	Matrix(int row, int col)
	{
		arr = new double[row * col];
		row_count = row;
		col_count = col;
	}
	double * el(int row, int col)
	{
		if (row >= row_count || col >= col_count) return 0;
		return arr + row + col * row_count;
	}
	~Matrix()
	{
		delete[] arr;
	}
	int row() {return row_count;}
	int col() {return col_count;}
	void operator*=(Matrix & A)
	{
		//Тут BLAS		(!!!)
		//			(!!!)
		//			(!!!)
		//			(!!!)
		//			(!!!)
		//			(!!!)
		//			(!!!)
		//Тут BLAS		(!!!)
		double * C = new double[row_count * col_count];
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_count, A.col_count, col_count, 1, arr, col_count, A.arr, A.col_count, 0, C, A.col_count);
		for (int i = 0; i < row_count * col_count; i++)
		{
			arr[i] = C[i];
		}
		delete[] C;
	}
};
int LUka(Matrix & L, Matrix & U)
{
	if (L.row() != U.row() || L.col() != U.col() || U.col() != L.row()) return 1;

	for (int i = 0; i < L.col(); i++)
	for (int k = 0; k < L.col(); k++)
		*L.el(i,k) = 0;

	for (int i = 0; i < U.col(); i++) *L.el(i,i) = 1;
	U *= L;
	
	for (int i = 0; i < U.col(); i++)
	{
		if (abs(*U.el(i,i)) < 0.0001) return 1;
		for (int k = i + 1; k < U.col(); k++)
		{
			double l = *U.el(k,i) / *U.el(i,i);
			*L.el(k,i) = l;
			for (int j = 0; j < U.col(); j++)
			{
				*U.el(k,j) -= *U.el(i,j) * l;
			}
		}
		if (abs(*U.el(i,i)) < 0.0001) return 1;
	}
	return 0;
}

istream & operator>>(istream & is, Matrix & m)
{
	for (int i = 0; i < m.row(); i++)
	{
		for (int k = 0; k < m.col(); k++)
		{
			is >> *m.el(i,k);
		}
	}
	return is;
}

ostream & operator<<(ostream & os, Matrix & m)
{
	os << endl;
	for (int i = 0; i < m.row(); i++)
	{
		for (int k = 0; k < m.col(); k++)
		{
			os << *m.el(i,k) << ' ';
		}
		os << endl;
	}
	return os;
}

int main()
{
	Matrix u(2,2), l(2,2);
	cin >> u;
    auto start = std::chrono::system_clock::now();
   	
	int err = LUka(l,u);

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
	if (err) cout << "it isn't possible" << endl;
	else cout << l << u << endl;
	
}

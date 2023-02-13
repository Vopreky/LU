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
	Matrix(Matrix & m)
	{
		int row = m.row_count;
		int col = m.col_count;
		arr = new double[row * col];
		row_count = row;
		col_count = col;
		for (int i = 0; i < row * col; i++)
			arr[i] = m.arr[i];
	}
	void load(Matrix & m)
	{
		for (int i = 0; i < row_count * col_count; i++)
			arr[i] = m.arr[i];
	}
	void clear()
	{
		for (int i = 0; i < row_count * col_count; i++)
			arr[i] = 0;
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
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_count, A.col_count, col_count, 1, A.arr, col_count, arr, A.col_count, 0, C, A.col_count);
		for (int i = 0; i < row_count * col_count; i++)
		{
			arr[i] = C[i];
		}
		delete[] C;
	}
	void operator-=(Matrix & A)
	{
		for (int i = 0; i < row_count * col_count; i++)
			arr[i] -= A.arr[i];
	}
};

int LUka(Matrix & L, Matrix & U, int l1, int l2, int u1, int u2)
{
//	cout << "LUKA_S" << endl;

	int b = L.row() / 2;
//	cout << b << endl;
//	cout << l1 << endl;
//	cout << l2 << endl;
//	cout << u1 << endl;
//	cout << u2 << endl;
//	cout << "LUKA_1" << endl;

	for (int i = 0; i < b; i++)
	for (int k = 0; k < b; k++)
		*L.el(i + b * l1, k + b * l2) = 0;
//	cout << "LUKA_2" << endl;

	for (int i = 0; i < b; i++) *L.el(i + b * l1, i + b * l2) = 1;
//	cout << "LUKA_3" << endl;
	
	for (int i = 0; i < b; i++)
	{
		//if (abs(*U.el(i + b * u1,i + b * u2)) < 0.0001) return 1;
		for (int k = i + 1; k < b; k++)
		{
			double l = *U.el(k + b * u1, i + b * u2) / *U.el(i + b * u1, i + b * u2);
			*L.el(k + b * l1, i + b * l2) = l;
			for (int j = 0; j < b; j++)
			{
				*U.el(k + b * u1, j + b * u2) -= *U.el(i + b * u1, j + b * u2) * l;
			}
		}
		//if (abs(*U.el(i + b * u1, i + b * u2)) < 0.0001) return 1;
	}
//	cout << "LUKA_E" << endl;
	return 0;
}

void LU(Matrix & L, Matrix & U)
{
	int b = L.row() / 2;
	U.clear();
	L.clear();
	//L.Get(0,1).clear();
	LUka(L,U,0,0,0,0);//LUka(L.Get(0,0), U.Get(0,0));
	for(int i = 0; i < b; i++)
	{
		for (int k = i + 1; k < b; k++)
		{
			double l = *L.el(k,i);//double l = *l11.el(k,i);
			for (int j = 0; j < b; j++)
			{
				*U.el(k,j + b) -= *U.el(i,j + b) * l;//*u12.el(k,j) -= *u12.el(i,j) * l;
			}
		}
	}
	//Matrix & l21 = L.Get(1,0);
	//Matrix & u11 = U.Get(0,0);
	//l21.load(U.Get(1,0));
	for (int i = 0; i < b; i++)
	{
		for (int j = 0; j < b; j++){
			*L.el(i + b,j) = *U.el(i + b,j);
		}
	}
	for (int i = 0; i < b; i++)
	{
		for (int k = i + 1; k < b; k++)
		{
			double l = *U.el(i,k) / *U.el(i,i);
			for (int j = 0; j < b; j++)
			{
				*L.el(j + b,k) -= *L.el(j + b,i) * l;//*l21.el(j,k) -= *l21.el(j,i) * l;
			}
		}

		for (int j = 0; j < b; j++)
		{
			*L.el(j + b,i) /= *U.el(i,i);//*l21.el(j,i) /= *u11.el(i,i);
		}
	}

	//U.Get(1,0).clear();
	//Matrix m(L.Get(1,0));
	Matrix m(b,b), m1(b,b);
	for (int i = 0; i < b; i++)
	{
		for (int j = 0; j < b; j++)
		{
			*m.el(i,j) = *L.el(i + b,j);
			*m1.el(i,j) = *U.el(i,j + b);
		}
	}
	//Вот тут использую blas, на этот раз с пользой			(!!!)
	//m *= U.Get(0,1);
	m *= m1;
	//Вот тут использую blas, на этот раз с пользой			(!!!)
	
	for (int i = 0; i < b; i++)
	{
		for (int j = 0; j < b; j++)
		{
			*U.el(i + b,j + b) -= *m.el(i,j);//U.Get(1,1) -= m;
		}
	}
	LUka(L,U,1,1,1,1);//LUka(L.Get(1,1), U.Get(1,1));
}

void test(Matrix & L, Matrix & U, Matrix & A)
{
	int sz = L.row();
	for (int i = 0; i < sz; i++)
	for (int k = 0; k < sz; k++)
	{
		double a = 0;
		for (int j = 0; j < sz; j++)
			a += (*L.el(i,j)) * (*U.el(j,k));
		if (abs(*(A.el(i,k)) - a) > 0.0001){
			cout << "error on (" << i << ", " << k << ")" << endl;
		}
	}
	for (int i = 0; i < sz; i++)
	for (int k = i + 1; k < sz; k++)
	{
		if (abs(*L.el(i,k)) > 0.0001){
			cout << "non-zero in L on (" << i << ", " << k << ")" << endl;
		}
		if (abs(*U.el(k,i)) > 0.0001){
			cout << "non-zero in L on (" << k << ", " << i << ")" << endl;
		}
	}
	for (int i = 0; i < sz; i++)
	{
		if (*L.el(i,i) != 1){
			cout << "non-one in L on (" << i << ", " << i << ")" << endl;
		}
	}
	cout << "OK" << endl;
}
istream & operator>>(istream & is, Matrix & M)
{
	for (int i = 0; i < M.row(); i++)
	{
		for (int k = 0; k < M.col(); k++)
		{
			is >> *M.el(i,k);
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
	const int b = 4;
	Matrix L(b,b), U(b,b), A(b,b);
	//cout << "Ввод матрицы блоками 11 12 21 22" << endl;
	cin >> U;
	/*
	for(int i = 0; i < b*2; i++)
		*U.el(i,i) = 1;
	for (int i = 0; i < 2; i++)
	for (int k = 0; k < 2; k++)
	A.Get(i,k).load(U.Get(i,k));
	*/
    auto start = std::chrono::system_clock::now();
	LU(L,U);
    auto end = std::chrono::system_clock::now();
	//cout << "L:\n" << L << endl;
	//cout << "U:\n" << U << endl;
   	


    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
    test(L,U,A);
//*/
}

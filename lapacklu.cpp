#include <iostream>
#include <chrono>
#include <ctime>
using namespace std;

extern "C"
{
extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
}
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
		double * C = new double[row_count * col_count];
		double alpha = 1;
		double beta = 0;
		
		//Тут теперь Лапак		(!!!)
		//				(!!!)
		//				(!!!)
		//				(!!!)
		//Тут теперь Лапак		(!!!)
		dgemm_("N", "N", &row_count, &A.col_count, &col_count, &alpha, arr, &col_count, A.arr, &A.col_count, &beta, C, &A.col_count);
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
class BigMrx
{
	Matrix A11;
	Matrix A12;
	Matrix A21;
	Matrix A22;
public:
	BigMrx(int s, int c): A11(s,c), A12(s,c), A21(s,c), A22(s,c) {}
	Matrix & Get(int i, int k)
	{
		if (i == 0 && k == 0)
			return A11;
		if (i == 0 && k == 1)
			return A12;
		if (i == 1 && k == 0)
			return A21;
		return A22;
	}
	double * el(int i, int k)
	{
		Matrix * m = 0;
		if (i >= A11.row())
		{
			i -= A11.row();
			if (k >= A11.col())
			{
				k -= A11.col();
				m = &A22;
			}
			else
			{
				m = &A21;
			}
		}
		else
		{
			if (k >= A11.col())
			{
				k -= A11.col();
				m = &A12;
			}
			else
			{
				m = &A11;
			}
		}
		return m->el(i,k);
	}
	int row()
	{
		return A11.row() * 2;
	}
	int col()
	{
		return A11.col() * 2;
	}
};

int LUka(Matrix & L, Matrix & U)
{
	if (L.row() != U.row() || L.col() != U.col() || U.col() != L.row()) return 1;

	for (int i = 0; i < L.col(); i++)
	for (int k = 0; k < L.col(); k++)
		*L.el(i,k) = 0;

	for (int i = 0; i < U.col(); i++) *L.el(i,i) = 1;
	
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

void LU(BigMrx & L, BigMrx & U)
{
	LUka(L.Get(0,0), U.Get(0,0));
	L.Get(0,1).clear();
	Matrix & u12 = U.Get(0,1);
	Matrix & l11 = L.Get(0,0);
	for(int i = 0; i < l11.col(); i++)
	{
		for (int k = i + 1; k < l11.row(); k++)
		{
			double l = *l11.el(k,i);
			for (int j = 0; j < u12.col(); j++)
			{
				*u12.el(k,j) -= *u12.el(i,j) * l;
			}
		}
	}
	Matrix & l21 = L.Get(1,0);
	Matrix & u11 = U.Get(0,0);
	l21.load(U.Get(1,0));
	for (int i = 0; i < u11.row(); i++)
	{
		for (int k = i + 1; k < u11.col(); k++)
		{
			double l = *u11.el(i,k) / *u11.el(i,i);
			for (int j = 0; j < l21.row(); j++)
			{
				*l21.el(j,k) -= *l21.el(j,i) * l;
			}
		}

		for (int j = 0; j < l21.row(); j++)
		{
			*l21.el(j,i) /= *u11.el(i,i);
		}
	}

	U.Get(1,0).clear();
	Matrix m(L.Get(1,0));
	
	//Эта функция теперь на лапаке держится
	m *= U.Get(0,1);
	//Эта функция теперь на лапаке держится
	
	U.Get(1,1) -= m;
	LUka(L.Get(1,1), U.Get(1,1));
}

void test(BigMrx & L, BigMrx & U, BigMrx & A)
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
	cout << "OK" << endl;
}
istream & operator>>(istream & is, BigMrx & M)
{
	for (int j1 = 0; j1 < 2; j1++)
	for (int j2 = 0; j2 < 2; j2++)
	{
		Matrix & m = M.Get(j1,j2);
		for (int i = 0; i < m.row(); i++)
		{
			for (int k = 0; k < m.col(); k++)
			{
				is >> *m.el(i,k);
			}
		}
	}
	return is;
}

ostream & operator<<(ostream & os, BigMrx & M)
{
	for (int j1 = 0; j1 < 2; j1++)
	for (int j2 = 0; j2 < 2; j2++)
	{
		os << endl;
		os << "Block(" << j1 << ", "<< j2 << ")" << endl;
		Matrix & m = M.Get(j1,j2);
		for (int i = 0; i < m.row(); i++)
		{
			for (int k = 0; k < m.col(); k++)
			{
				os << *m.el(i,k) << ' ';
			}
			os << endl;
		}
	}
	return os;
}

int main()
{
	BigMrx L(2,2), U(2,2), A(2,2);
	cin >> U;
	for (int i = 0; i < 2; i++)
	for (int k = 0; k < 2; k++)
	A.Get(i,k).load(U.Get(i,k));
    auto start = std::chrono::system_clock::now();
	LU(L,U);
    auto end = std::chrono::system_clock::now();
	cout << "L:\n" << L << endl;
	cout << "U:\n" << U << endl;
   	


    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
    test(L,U,A);
//*/
}

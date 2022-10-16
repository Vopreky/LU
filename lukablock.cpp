#include <iostream>
#include <chrono>
#include <ctime>
using namespace std;


const int b = 32;
void init_c(double * block)
{
	for (int i = 0; i < b; i++)
	for (int k = 0; k < b; k++)
		block[i*b + k] = 0;
}

void init_a(double * arr, int x, int y, int sz, double * block)
{
	for (int i = 0; i < b; i++)
	for (int k = 0; k < b; k++)
		block[i*b + k] = arr[(x + i) * sz + (y + k)];
}

void small_mult(double * A, double * B, double * res)
{
/*
	cout << "mult small matrixes" << endl;
	print(A,b);
	cout << endl;
	print(B,b);
	cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
	cout << endl;
	cout << endl;
*/
	for (int i = 0; i < b; i++)
	for (int k = 0; k < b; k++)
	{
		res[i*b + k] = 0;
		for (int j = 0; j < b; j++)
		{
			res[i*b + k] += A[j*b + k]*B[i*b + j];
		}
	}
}

void load_block(double * A, double * block, int x, int y, int sz)
{
	for (int i = 0; i < b; i++)
	for (int k = 0; k < b; k++)
	{
		A[(x + i)*sz + (y + k)] = block[i * b + k];
	}
}


class MRX
{

	double * arr;
	int sz;
public:
	MRX(int sz, double val = 0)
	{
		if (sz % b != 0)
			throw "Bad sz";
		this->sz = sz;
		arr = new double[sz*sz];
		for (int i = 0; i < sz*sz; i++) arr[i] = val;
	}
	double * GetBlock(int s, int c)
	{
		int n = sz / b;
		int offset = (c * n + s) * b*b;
		return arr + offset;
	}
	double * el(int s, int c)
	{
		int block_s = s / b;
		int block_c = c / b;
		double * block = GetBlock(block_s, block_c);
		c -= block_c * b;
		s -= block_s * b;
		int offset = (c * b + s);
		return block + offset;
	}
	void print()
	{
		for (int i = 0; i < sz; i++)
		{
			for (int k = 0; k < sz; k++)
			{
				cout << *el(i,k) << ' ';
			}
			cout << endl;
		}
	}
	int size()
	{
		return sz;
	}
~MRX()
{
	delete[] arr;
}
};
MRX * mult(MRX & A, MRX & B, int sa, int ca, int sb, int cb, int sz)
{
	MRX * res = new MRX(sz,0);
	//return res;

	double * block_a;
	double * block_b;
	double * block_c;
	
	for (int i = 0; i < sz; i+=b)
	for (int j = 0; j < sz; j+=b)
	{
		block_c = res->GetBlock(i/b,j/b);
		init_c(block_c);
		for (int k = 0; k < sz; k+=b)
		{
			block_a = A.GetBlock(i/b + sa/b, k/b + ca/b);
			block_b = B.GetBlock(k/b + sb/b, j/b + cb/b);
			small_mult(block_a, block_b, block_c);
		}
	}
	return res;
}


void rec(MRX & A, MRX & L, MRX & U, int sa, int ca, int sl, int cl, int su, int cu, int sz)
{
	if (sz == b)
	{
		double * u_block = U.GetBlock(su/b,cu/b);
		double * a_block = A.GetBlock(sa/b,ca/b);
		for (int i = 0 ; i < b*b; i++)
			u_block[i] = a_block[i];
		for (int i = 0; i < b; i++) *L.el(sl + i, cl + i) = 1;
        
		for (int i = 0; i < b; i++)
		{
			if (abs(*U.el(su + i, cu + i)) < 0.0001) throw "error minor ~ 0";
			for (int k = i + 1; k < b; k++)
			{
				double l = *U.el(su + k, cu + i) / *U.el(su + i, cu + i);
				*L.el(sl + k, cl + i) = l;
				for (int j = 0; j < b; j++)
				{
					*U.el(su + k, cu + j) -= *U.el(su + i, cu + j) * l;
				}
			}
		}
	}
	else if (sz > b)
	{
		rec(A, L, U, sa, ca, sl, cl, su, cu, sz / 2);
	//error is here	
		for (int i = 0; i < sz/2 - 1; i++)
		{
			for (int k = i + 1; k < sz/2; k++)
			{
				double l = (*L.el(sl + k, cl + i));
				for (int j = 0; j < sz/2; j++)
				{
					*A.el(sa + k, ca + sz/2 + j) -= *A.el(sa + i, ca + sz/2 + j) * l;
				}
			}
		}

		for (int i = 0; i < sz/2/b; i++)
		for (int k = 0; k < sz/2/b; k++)
		{
			double * block_u = U.GetBlock(su/b + i, cu/b + sz/2/b + k);
			double * block_a = A.GetBlock(su/b + i, cu/b + sz/2/b + k);
			for (int j = 0; j < b*b; j++)
				block_u[j] = block_a[j];
		}

		for (int i = 0; i < sz/2 - 1; i++)
		{
			for (int k = i + 1; k < sz/2; k++)
			{
				double l = (*U.el(su + i, cu + k) / *U.el(su + i, cu + i));
				for (int j = 0; j < sz/2; j++)
				{
					*A.el(sa + sz/2 + j, ca + k) -= *A.el(sa + sz/2 + j, ca + i) * l;
				}
			}
		}
		
		for (int i = 0; i < sz/2; i++)
		{
			for (int k = 0; k < sz/2; k++)
			{
				*A.el(sa + sz/2 + k, ca + i) /= *U.el(su + i, cu + i);
			}
		}
		
		for (int i = 0; i < sz/2/b; i++)
		for (int k = 0; k < sz/2/b; k++)
		{
			double * block_l = L.GetBlock(su/b + sz/2/b + i, cu/b + k);
			double * block_a = A.GetBlock(su/b + sz/2/b + i, cu/b + k);
			for (int j = 0; j < b*b; j++)
				block_l[j] = block_a[j];
		}
	//error is here	

		MRX * l21u12 = mult(L,U, sl + sz / 2, cl, su, cu + sz / 2, sz / 2); 
		//cout << "l21u12 = (" << sl << ", " << sz/2 << ")" << endl;
		//l21u12->print();
		//cout << endl;
		//L.print();
		//cout << endl;
		//U.print();
		//cout << endl;
		for (int i = 0; i < sz / 2 / b; i++)
		for (int k = 0; k < sz / 2 / b; k++)
		{
			double * block_aa = A.GetBlock(sa/b + sz/2/b + i, ca/b + sz/2/b + k);
			double * block_lu = l21u12->GetBlock(i,k);
			for (int j = 0; j <= b*b; j++)
			{
				//cout << "aa " << block_aa[j] << endl;
				//cout << "lu " << block_lu[j] << endl;
				block_aa[j] -= block_lu[j];
			}
		}
		delete l21u12;
		rec(A, L, U, sa + sz/2, ca + sz/2, sl + sz/2, cl + sz/2, su + sz/2, cu + sz/2, sz/2);
	}
	else throw "What???";
}

void LUkablock(MRX & A, MRX & L, MRX & U)
{
	rec(A,L,U,0,0,0,0,0,0,A.size());
}

int main(int argc, char ** argv)
{
	MRX m1(128*b,1), m2(128*b), m3(128*b);

	for (int i = 0; i < 128*b; i++)
		*m1.el(i,i) = 2;

    auto start = std::chrono::system_clock::now();
   	
	try
	{
		LUkablock(m1,m2,m3);
	}
	catch (const char * s)
	{
		cout << s << endl;
		return 0;
	}

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
	
}




























#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
using namespace std;

extern "C"
{
extern void dgemm_(const char* transa,
		const char * transb,
		int * m,
		int * n,
		int * k,
		double * alpha,
		double * A,
		int * lda,
		double * B,
		int * ldb,
		double * beta,
		double * C,
		int * ldc);
extern void dstevd_(const char * jobz,
		int * n,
		double * d,
		double * e,
		double * z,
		int * ldz,
		double * work,
		int * lwork,
		int * iwork,
		int * liwork,
		int * info);
}

double norma(double * v, int n){
	double a = 0;
	for (int i = 0; i < n; i++){
		a += v[i]*v[i];
	}
	return sqrt(a);
}

inline void sum(double * a, double * b, double * c, int n){
	for (int i = 0; i < n; i++)
		c[i] = a[i] + b[i];
}
inline void mult(double * a, double b, double * c, int n){
	for (int i = 0; i < n; i++)
		c[i] = a[i] * b;
}
inline void load(double * a, double * c, int n){
	for (int i = 0; i < n; i++)
		c[i] = a[i];
}
inline void nullify(double * a, int n){
	for (int i = 0; i < n; i++)
		a[i] = 0;
}
//A - matrix we work with
//n - sz of A dims
//matvec -> if (revFlag) {y = Ax y - (n, 1) x - (m, 1)} 
//		    else {y = xA y - (1, m) x - (1, n)} 
//		    A - (n, m)
void lanc(double * A, int n, void(*matvec) (double * A, double * x, double * y, int n, int m, bool revFlag), int maxK = 10){
	double * workArray = new double[n];
	double b;
	vector<double*> sys;//vector<mrx> sys;
	vector<double> aArray;
	vector<double> bArray;
	double * q = new double[n];//mrx q(A.s);
	nullify(q,n);

	q[0] = 1;//q.set(1);
	sys.push_back(q);
	q = new double[n];
	load(sys[0], q, n);
	matvec(A, q, workArray, n, n, true);
	matvec(workArray, q, workArray, 1, n, false);
	double a = -workArray[0];//double a = -((q.T() * A * q).get());
	aArray.push_back(a);
	matvec(A, q, workArray, n, n, false);
	mult(q, a, q, n);
	sum(workArray, q, q, n);
	//q = A * q + q * a;
	int k = 1;
	while (norma(q, n) > 0.001 && k < maxK){
		int rank =  omp_get_thread_num();
		mult(q, 1 / norma(q, n), q, n);
		//q = q / norma(q);
		sys.push_back(q);
		q = new double[n];
		k++;
		
		matvec(A, sys[k-1], workArray, n, n, true);
		matvec(workArray, sys[k-1], workArray, 1, n, false);	
		a = -workArray[0];//-((sys[k-1].T() * A * sys[k-1]).get());
		matvec(A, sys[k-1], workArray, n, n, true);
		matvec(workArray, sys[k-2], workArray, 1, n, false);
		b = -workArray[0];//-((sys[k-1].T() * A * sys[k-2]).get());
		
		aArray.push_back(a);
		bArray.push_back(b);

		matvec(A, sys[k-1], q, n, n, false);
		
		
		#pragma omp parallel 
		{
			double * my_q = new double [n];
			for (int i = 0; i < n; i++)
				my_q[i] = 0;
			#pragma omp for
			for (int i = 0; i < n; i++){
				my_q[i] = sys[k-1][i] * a + sys[k-2][i] * b;
			}
			#pragma omp critical
			for (int i = 0; i < n; i++){
				q[i] += my_q[i];
			}
			delete[] my_q;
			//mult(sys[k-1], a, workArray, n);
			//sum(q, workArray, q, n);
			//mult(sys[k-2], b, workArray, n);
			//sum(q, workArray, q, n);
		}
		//q = A * sys[k-1] + sys[k-1] * a + sys[k-2] * b;
	}
	cout << "k = " << k << endl;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < k; j++){
			cout << sys[j][i] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for (auto p = aArray.begin(); p != aArray.end(); p++){
		cout << *p << " ";
	}
	cout << endl;
	for (auto p = bArray.begin(); p != bArray.end(); p++){
		cout << *p << " ";
	}
	cout << endl;
	double * mid = new double[k];
	double * bot = new double[k - 1];
	for (int i = 0; i < k; i++){
		mid[i] = -aArray[i];
	}
	for (int i = 0; i < k - 1; i++){
		bot[i] = -bArray[i];
	}
	double * arr1 = new double[k];
	int * arr2 = new int[k];
	int k2 = 4*k;
	int info;
	int one = 1;
	dstevd_("N", &k, mid, bot, 0, &one, arr1, &k, arr2, &k, &info);
	//dgeev_("N", "N", &k, T, &k, arr1, arr2, arr3, &k, arr4, &k, arr5, &k2, &info);
	cout << "info = " << info << endl;
	cout << "result:" << endl;
	for (int i = 0; i < k; i++){
		cout << mid[i] << endl;
	}

delete[] arr1;
delete[] arr2;
delete[] q;
delete[] workArray;
delete[] mid;
delete[] bot;
for (auto i = sys.begin(); i != sys.end(); i++)
	delete[] *i;
}

void matvec(double * A, double * x, double * y, int n, int m, bool revFlag){
	double alpha = 1, beta = 0;
	int one = 1;
	if (!revFlag){
		double * temp = new double[n];
		dgemm_("N", "N", &n, &one, &m, &alpha, A, &n, x, &m, &beta, temp, &n);
		load(temp, y, n);
		delete[] temp;
	}
	else{
		double * temp = new double[m];
		dgemm_("N", "N", &one, &m, &n, &alpha, x, &one, A, &n, &beta, temp, &one);
		load(temp, y, m);
		delete[] temp;
	}
}

int main(){
	int n = 25000;
	double * A = new double[n*n];
	for (int i = 0; i < n*n; i++)
		A[i] = 1;
	lanc(A,n,matvec);
	delete[] A;
}



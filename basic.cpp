#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;
extern "C"
{
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

template<class T>
void vecLoad(T * vec, T * res, int size){
	for (int i = 0; i < size; i++)
		res[i] = vec[i];
}
template<class T>
void vecMult(T * vec, T a, T * res, int size){
	for (int i = 0; i < size; i++)
		res[i] = vec[i] * a;
}
template<class T>
void vecSum(T * v1, T * v2, T * res, int size){
	for (int i = 0; i < size; i++)
		res[i] = v1[i] + v2[i];
}
template<class T>
T dotprod(T * v1, T * v2, int size){
	T res = 0;
	for (int i = 0; i < size; i++){
		res += v1[i] * v2[i];
	}
	return res;
}
template<class T>
void matvec(T * mat, T * vec, T * res, int size){
	for (int i = 0; i < size; i++){
		res[i] = 0;
		for (int j = 0; j < size; j++){
			res[i] += mat[i * size + j] * vec[j];
		}
	}
}
template<class T>
T norma(T * v, int size){
	return sqrt(dotprod(v, v, size));
}
template<class T>
T * lanc(T * A, T * b, int k, int & size){
	int count = 1;
	vector<T*> q;
	vector<T> a;
	vector<T> g;

	T * ptr = new T[size];
	T * v1 = new T[size];
	T * v2 = new T[size];
	vecMult(ptr, (T)0, ptr, size);
	q.push_back(ptr);
	ptr = new T[size];
	vecMult(b, 1 / norma(b, size), ptr, size);
	q.push_back(ptr);

	a.push_back(0);
	g.push_back(0);
	for (int j = 1; j <= k; j++){
		T * z = new T[size];
		matvec(A, q[j], z, size);
		a.push_back(dotprod(q[j], z, size));
		count++;
		vecMult(q[j-1], -g[j-1], v1, size);
		vecMult(q[j], -a[j], v2, size);
		vecSum(z, v1, z, size);
		vecSum(z, v2, z, size);
		for (int i = 1; i < j; i++){
			vecMult(q[i], -dotprod(z, q[j], size), v1, size);
			vecSum(z, v1, z, size);
		}
		for (int i = 1; i < j; i++){
			vecMult(q[i], -dotprod(z, q[j], size), v1, size);
			vecSum(z, v1, z, size);
		}
		g.push_back(norma(z, size));
		if (g[j] == 0)
			break;
		vecMult(z, 1 / g[j], z, size);
		q.push_back(z);
	}
	size = count;
	T * res = new T[count * count];
	vecMult(res, (T)0, res, size * size);
	for (int i = 0; i < count; i++){
		res[count * i + i] = a[i + 1];
	}
	for (int i = 0; i < count - 1; i++){
		res[count * i + (i + 1)] = g[i + 2];
		res[count * (i + 1) + i] = g[i + 2];
	}
	for (auto i = q.begin(); i != q.end(); i++){
		delete[] *i;
	}
	delete[] v1;
	delete[] v2;
	return res;
}


int main(int argc, char ** argv){
	const int size = 12025;
	double * A = new double[size * size];
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			A[i * size + j] = 1;
	
	double b[size];
	for (int i = 0; i < size; i++)
		b[i] = i + 1;
	int k = 25;

	int mrxSize = size;
	auto start = high_resolution_clock::now();
	double * T = lanc(A,b,k,mrxSize);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;
	
	
	double * mid = new double[mrxSize];
	double * bot = new double[mrxSize - 1];
	for (int i = 0; i < mrxSize; i++){
		mid[i] = T[mrxSize * i + i];
	}
	for (int i = 0; i < mrxSize - 1; i++){
		bot[i] = T[mrxSize * i + i + 1];
	}
	double * arr1 = new double[mrxSize];
	int * arr2 = new int[mrxSize];
	int k2 = 4*mrxSize;
	int info;
	int one = 1;
	dstevd_("N", &mrxSize, mid, bot, 0, &one, arr1, &mrxSize, arr2, &mrxSize, &info);
	cout << "info = " << info << endl;
	cout << endl << "самое большое сз матрицы = " << size << ", результате должен быть что-то близкое" << endl << endl;
	cout << "result:" << endl;
	cout << mid[mrxSize - 2] << endl;
	delete[] mid;
	delete[] bot;
	delete[] T;
	delete[] A;
}

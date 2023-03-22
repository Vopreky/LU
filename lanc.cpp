#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

extern "C"
{
extern void dgeev_(const char*, const char*, int, double *, int, double *, double *, double *, int, double *, int, int *);
}

class mrx{
public:
	double * arr;
	int s,c;
	mrx(int s, int c = 1){
		arr = new double[s*c];
		this->s = s;
		this->c = c;
	}
	mrx(const mrx & m){
		c = m.c;
		s = m.s;
		arr = new double[s*c];
		for (int i = 0; i < s * c; i++)
			arr[i] = m.arr[i];
	}
	const mrx operator*(double a) const{
		mrx res(s,c);
		for (int i = 0; i < s*c; i++)
			res.arr[i] = arr[i] * a;
		return res;
	}
	const mrx operator/(double a) const{
		return (*this) * (1/a);
	}
	const mrx & operator=(const mrx & m){
		delete[] arr;
		c = m.c;
		s = m.s;
		arr = new double[s*c];
		for (int i = 0; i < s * c; i++)
			arr[i] = m.arr[i];
		return *this;
	}
	const mrx operator*(const mrx & m) const{
		if (c != m.s)
			throw "ET1";
		mrx res(s,m.c);
		for (int i = 0; i < s; i++){
			for (int j = 0; j < m.c; j++){
				for (int k = 0; k < c; k++){
					res.arr[i * res.c + j] += arr[i * c + k] * m.arr[k * m.c + j];
				}
			}
		}
		return res;
	}
	const mrx operator+(const mrx & m) const{
		if (c != m.c)
			throw "ET2";
		if (s != m.s)
			throw "ET3";
		mrx res(s,c);
		for (int i = 0; i < s * c; i++)
			res.arr[i] = arr[i] + m.arr[i];
		return res;
	}
	const mrx operator-(const mrx & m) const{
		if (c != m.c)
			throw "ET2";
		if (s != m.s)
			throw "ET3";
		mrx res(s,c);
		for (int i = 0; i < s * c; i++)
			res.arr[i] = arr[i] - m.arr[i];
		return res;
	}
	~mrx(){
		delete[] arr;
	}
	void print() const{
		for (int i = 0; i < s; i++){
			for (int j = 0; j < c; j++){
				cout << arr[i * c + j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	double get(int i = 0, int j = 0) const{
		return arr[i * c + j];
	}
	void set(double a, int i = 0, int j = 0){
		arr[i * c + j] = a;
	}
	const mrx T() const{
		mrx res(c,s);
		for (int i = 0; i < s; i++)
		for (int j = 0; j < c; j++)
			res.set(get(i,j),j,i);
		return res;
	}
};
double norma(const mrx & v){
	mrx a = v.T() * v;
	return sqrt(a.get());
}

void lanc(const mrx & A, int maxK = 10){
	double b;
	vector<mrx> sys;
	vector<double> aArray;
	vector<double> bArray;
	mrx q(A.s);

	q.set(1);
	sys.push_back(q);
	double a = -((q.T() * A * q).get());
	aArray.push_back(a);
	q = A * q + q * a;
	int k = 1;
	while (norma(q) > 0.001 && k < maxK){
		q = q / norma(q);
		sys.push_back(q);
		k++;
		a = -((sys[k-1].T() * A * sys[k-1]).get());
		b = -((sys[k-1].T() * A * sys[k-2]).get());
		aArray.push_back(a);
		bArray.push_back(b);
		q = A * sys[k-1] + sys[k-1] * a + sys[k-2] * b;
	}
	cout << "k = " << k << endl;
	for (int i = 0; i < A.s; i++){
		for (int j = 0; j < k; j++){
			cout << sys[j].get(i) << " ";
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
	mrx T(k,k);
	for (int i = 0; i < k; i++){
		T.set(aArray[i],i,i);
	}
	for (int i = 0; i < k - 1; i++){
		T.set(bArray[i],i + 1, i);
		T.set(bArray[i],i, i + 1);
	}
	T.print();
	double * arr1 = new double[k*k];
	double * arr2 = new double[k*k];
	double * arr3 = new double[k*k];
	double * arr4 = new double[k*k];
	int info;
	dgeev_("N", "N", k, T.arr, k, arr1, arr2, arr3, 0, arr4, 0, &info);
	cout << "info = " << info << endl;
	for (int i = 0; i < k; i++){
		cout << arr1[i] << " + i * " << arr2[i] << endl;
	}

delete[] arr1;
delete[] arr2;
delete[] arr3;
delete[] arr4;
}

int main(){
	mrx A(3,3);
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		A.set(1, i, j);
	lanc(A);
	/*
	mrx B(3,1);
	for (int i = 0; i < 3; i++)
		B.set(1,i);
	mrx C = A * B;
	A.print();
	B.print();
	C.print();
	C.T().print();
	(C + B).print();
	*/
}



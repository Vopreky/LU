program main
implicit none (type, external)
external :: dgemm
external :: dtrsm
double precision   :: L11(2, 2)
double precision   :: L12(2, 2)
double precision   :: L21(2, 2)
double precision   :: L22(2, 2)
       
double precision   :: U11(2, 2)
double precision   :: U12(2, 2)
double precision   :: U21(2, 2)
double precision   :: U22(2, 2)
      
double precision   :: A11(2, 2)
double precision   :: A12(2, 2)
double precision   :: A21(2, 2)
double precision   :: A22(2, 2)

double precision :: one;
double precision :: two;
        
double precision:: L(2,2)
double precision:: U(2,2)
integer :: i, j, k, n
real :: ll
one = 1;
two = 2;
! Шаг 1 - ввод матрицы 

A11(1,1) = 1
A11(1,2) = 2
A11(2,1) = 3
A11(2,2) = 4

A12(1,1) = 1
A12(1,2) = 0
A12(2,1) = 0
A12(2,2) = 1

A21(1,1) = 1
A21(1,2) = 0
A21(2,1) = 0
A21(2,2) = 1

A22(1,1) = 1
A22(1,2) = 2
A22(2,1) = 3
A22(2,2) = 4

U11 = A11;
U12 = A12;
U21 = A21;
U22 = A22;

! Шаг 2 - вычисление L11, U11

L = L11;
U = U11;
n = 2;
        L = L * 0;
        i = 1;
        do while (i <= n)
                L(i, i) = 1;
                i = i + 1;
        end do
        i = 1;
        do while(i <= n)
                k = i + 1;
                do while(k <= n)
                        ll = U(k, i) / U(i, i);
                        L(k, i) = ll;
                        j = 1;
                        do while(j <= n)
                                U(k, j) = U(k, j) - U(i, j) * ll;
                        j = j + 1
                        end do
                k = k + 1;
                end do
        i = i + 1;
        end do
L11 = L;
U11 = U;

! Шаг 3 - вычисление L21, U12

!call dtrsm(side, uplo, trA, diag, m, n, alpha, A, LDA, B, LDB)
U12 = A12;
L21 = A21;
call dtrsm('L', 'L', 'N', 'N', 2, 2, one, L11, 2, U12, 2);
call dtrsm('R', 'U', 'N', 'N', 2, 2, one, U11, 2, L21, 2);

! Шаг 4 - вычисление L22, U22

one = -1;
two = 1;
        print *,U22
        print *,L21
        print *,''
call dgemm('N', 'N', 2, 2, 2, one, L21, 2, U12, 2, two, U22, 2)
        print *,U22
        print *,L21
        print *,''
L = L22;
U = U22;
n = 2;
L = L * 0;
i = 0;
do while (i < n)
        L(i + 1,i + 1) = 1;
        i = i + 1;
end do
i = 0;
do while(i < n)
        k = i + 1;
        do while(k < n)
                ll = U(k + 1, i + 1) / U(i + 1, i + 1);
                L(k + 1, i + 1) = ll;
                j = 0;
                do while(j < n)
                        U(k + 1, j + 1) = U(k + 1, j + 1) - U(i + 1, j + 1) * ll;
                j = j + 1
                end do
        k = k + 1;
        end do
i = i + 1;
end do
L22 = L;
U22 = U;

L12 = L12 * 0;
u21 = u21 * 0;

print *, A11
print *, A12
print *, A21
print *, A22
print *,' '
print *, L11
print *, L12
print *, L21
print *, L22
print *,' '
print *, U11
print *, U12
print *, U21
print *, U22

end program main

!subroutine dgemm	(	character 	TRANSA,
!character 	TRANSB,
!integer 	M,
!integer 	N,
!integer 	K,
!double precision 	ALPHA,
!double precision, dimension(lda,*) 	A,
!integer 	LDA,
!double precision, dimension(ldb,*) 	B,
!integer 	LDB,
!double precision 	BETA,
!double precision, dimension(ldc,*) 	C,
!integer 	LDC
!)

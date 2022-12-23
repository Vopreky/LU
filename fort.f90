SUBROUTINE luka(L, U, n)
        real L(0:3);
        real U(0:3);
        integer  i, j, k, n;
        real ll
        L = L * 0;
        i = 0;
        do while (i < n)
                L(i * n + i) = 1;
                i = i + 1;
        end do
        i = 0;
        do while(i < n)
                k = i + 1;
                do while(k < n)
                        ll = U(k * n + i) / U(i * n + i);
                        L(k * n + i) = ll;
                        j = 0;
                        do while(j < n)
                                U(k * n + j) = U(k * n + j) - U(i * n + j) * ll;
                        j = j + 1
                        end do
                k = k + 1;
                end do
        i = i + 1;
        end do
!int LUka(Matrix & L, Matrix & U)
!{
!	if (L.row() != U.row() || L.col() != U.col() || U.col() != L.row()) return 1;
!
!	for (int i = 0; i < U.col(); i++) *L.el(i,i) = 1;
!
!	for (int i = 0; i < U.col(); i++)
!	{
!		if (abs(*U.el(i,i)) < 0.0001) return 1;
!		for (int k = i + 1; k < U.col(); k++)
!		{
!			double l = *U.el(k,i) / *U.el(i,i);
!			*L.el(k,i) = l;
!			for (int j = 0; j < U.col(); j++)
!			{
!				*U.el(k,j) -= *U.el(i,j) * l;
!			}
!		}
!		if (abs(*U.el(i,i)) < 0.0001) return 1;
!	}
!	return 0;
!}

END SUBROUTINE luka
IMPLICIT none
REAL d
REAL A(0:3);
REAL B(0:3);
REAL C(0:3);
INTEGER n,i,k,j;
REAL val;
val = 0;
n = 2;
B(0) = 1;
B(1) = 2;
B(2) = 3;
B(3) = 4;
C = B;
CALL luka(A,B,n)
print *, A
print *, B
i = 0;
k = 0;
do while (i < n)
        do while (k < n)
                j = 0;
                d = 0;
                do while (j < n)
                        d = d + A(i*n + j) * B(j*n + k);
                j = j + 1;
                end do
                if (abs(d - C(i*n + k)) > 0.001) then
                        print *,"error i guess";
                end if
        k = k + 1;
        end do
i = i + 1;
end do
print *,"Ok?"
end

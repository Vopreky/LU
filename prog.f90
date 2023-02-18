
subroutine solve1(L,U,A)
        integer, parameter :: sz = 2
        external :: dtrsm
        integer :: i = 0
        integer :: j = 0
        double precision :: L(1:sz * 2, 1:sz * 2)
        double precision :: U(1:sz * 2, 1:sz * 2)
        double precision :: A(1:sz * 2, 1:sz * 2)
        double precision :: L11(1:sz, 1:sz)
        double precision :: U12(1:sz, 1:sz)
        double precision :: one = 1

        do i = 1,sz
                do j = 1,sz
                        U12(i,j) = A(i,j + sz)
                        L11(i,j) = L(i,j)
                end do
        end do
        call dtrsm('L', 'L', 'N', 'N', sz, sz, one, L11, sz, U12, sz);
        do i = 1,sz
                do j = 1,sz
                        U(i,j + sz) = U12(i,j)
                end do
        end do
        
        end subroutine

subroutine solve2(L,U,A)
        integer, parameter :: sz = 2
        external :: dtrsm
        integer :: i = 0
        integer :: j = 0
        double precision :: L(1:sz * 2, 1:sz * 2)
        double precision :: U(1:sz * 2, 1:sz * 2)
        double precision :: A(1:sz * 2, 1:sz * 2)
        double precision :: L21(1:sz, 1:sz)
        double precision :: U11(1:sz, 1:sz)
        double precision :: one = 1

        do i = 1,sz
                do j = 1,sz
                        U11(i,j) = U(i,j)
                        L21(i,j) = A(i + sz,j)
                end do
        end do
        call dtrsm('R', 'U', 'N', 'N', sz, sz, one, U11, sz, L21, sz);
        do i = 1,sz
                do j = 1,sz
                        L(i + sz,j) = L21(i,j)
                end do
        end do
        
        end subroutine

subroutine LU(L,U,A)
        integer, parameter :: sz = 2
        integer :: i = 0
        integer :: j = 0
        integer :: k = 0
        double precision :: L(1:sz * 2, 1:sz * 2)
        double precision :: U(1:sz * 2, 1:sz * 2)
        double precision :: A(1:sz * 2, 1:sz * 2)
        double precision :: T(1:sz, 1:sz)
        double precision :: x = 0

        L = 0
        U = 0
        do i = 1, sz
        L(i,i) = 1;
                do j = i, sz   
                        x = 0  
                        do k = 1,i-1
                                x = x + L(i,k) * U(k,j)
                        end do 
                        U(i,j) = A(i,j) - x
                end do   
                
                do j = i + 1, sz     
                        x = 0;
                        do k = 1,i-1
                                x = x + L(j,k) * U(k,i)
                        end do
                        L(j,i) = (A(j,i) - x)/U(i,i)
                end do   
        end do

        call solve1(L,U,A)
        call solve2(L,U,A)

        do i = 1,sz
                do j = 1,sz
                        x = 0
                        do k = 1,sz
                                x = x + L(sz + i,k)*U(k,j + sz)
                        end do
                        A(i + sz, j + sz) = A(i + sz, j + sz) - x
                end do
        end do
        do i = 1, sz
        L(i + sz,i + sz) = 1;
                do j = i, sz   
                        x = 0  
                        do k = 1,i-1
                                x = x + L(i + sz,k + sz) * U(k + sz,j + sz)
                        end do 
                        U(i + sz,j + sz) = A(i + sz,j + sz) - x
                end do   
                
                do j = i + 1, sz
                        x = 0;
                        do k = 1,i-1
                                x = x + L(j + sz, k + sz) * U(k + sz,i + sz)
                        end do
                        L(j + sz,i + sz) = (A(j + sz,i + sz) - x)/U(i + sz,i + sz)
                end do   
        end do

        end subroutine
program main
        integer, parameter :: sz = 2
        double precision :: L(1:sz * 2, 1:sz * 2)
        double precision :: U(1:sz * 2, 1:sz * 2)
        double precision :: A(1:sz * 2, 1:sz * 2)
        integer :: row = 0
        integer :: col = 0
        A = 1
        do row = 1,sz*2
                A(row,row) = 2;
        end do
        do row=1, sz*2
        write(*,*) (A(row,col),col=1,sz*2)
        enddo
        print *, ""
        call LU(L,U,A)

        do row=1, sz*2
        write(*,*) (L(row,col),col=1,sz*2)
        enddo
        print *, ""
        do row=1, sz*2
        write(*,*) (U(row,col),col=1,sz*2)
        enddo
        end program

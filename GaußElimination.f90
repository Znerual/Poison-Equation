    MODULE LinearSystems
    IMPLICIT NONE
    private
    public lufaktorisierungGauss
    public solveGauss
    public solveThomas
    public showMatrix
    INTEGER :: i,j,k,maxi,l,ii,jj
    DOUBLE PRECISION :: start_time, stop_time

    contains
    subroutine lufaktorisierungGauss(matrix, l_matrix,bb, n, m, lutime )
    DOUBLE PRECISION, DIMENSION(:,:), intent(inout)  :: matrix, l_matrix
    DOUBLE PRECISION, DIMENSION(:), intent(inout)    :: bb
    INTEGER, intent(in) :: n,m
    double precision,intent(out) :: lutime

    DOUBLE PRECISION, DIMENSION(:),allocatable    :: tmpmatrix1,tmpmatrix2
    DOUBLE PRECISION, DIMENSION(:,:), allocatable  :: ll

    CALL cpu_time(start_time)

    allocate(tmpmatrix1(m))
    allocate(tmpmatrix2(m))
    allocate(ll(n,m))
    tmpmatrix1(:)=0.d0
    tmpmatrix2(:)=0.d0
    l_matrix(:,:)=0.d0
    ll(:,:)=0.d0
    ! reduce generic matrix M to the upper triangular matrix U
    l=1
    k=1
    DO ii=1,n
        DO jj=k,m
            j=jj
            maxi=j
            DO i=j+1,n
                IF(ABS(matrix(i,j)).GT.matrix(maxi,j))THEN
                    maxi=i
                END IF
            END DO
            IF(matrix(k,k).NE.0)THEN
                DO i=l+1,n
                    ll(i,k)=matrix(i,k)/matrix(k,k)
                    DO j=k,n
                        matrix(i,j)=matrix(i,j)-ll(i,k)*matrix(k,j)
                    END DO
                    bb(i)=bb(i)-ll(i,k)*bb(k)
                END DO
                l=l+1
                k=k+1
            END IF
        END DO
    END DO
    ! store the lower triangular matrix L
    DO i=1,n
        DO j=1,i
            IF(i.EQ.j)l_matrix(i,j)=1.d0
            IF(i.NE.j)l_matrix(i,j)=ll(i,j)
        END DO
    END DO
    CALL cpu_time(stop_time)
    lutime = stop_time - start_time

    deallocate(tmpmatrix1)
    deallocate(tmpmatrix2)
    deallocate(ll)
    end subroutine lufaktorisierungGauss

    subroutine solveGauss(bb, matrix, n,m, solvetime)
    DOUBLE PRECISION, DIMENSION(:), intent(inout)    :: bb
    DOUBLE PRECISION, DIMENSION(:,:), intent(in)  :: matrix
    INTEGER, intent(in) :: n,m
    DOUBLE PRECISION, intent(out) :: solvetime

    CALL cpu_time(start_time)
    ! solve by backward substitution Ux=z (here x is overwritten on b and bb is z)
    DO ii=1,m
        i=m+1-ii
        DO j=i+1,n
            bb(i) = bb(i)-matrix(i,j)*bb(j)
        END DO
        IF(matrix(i,i).NE.0)THEN
            bb(i)=bb(i)/matrix(i,i)
        ELSE
            WRITE(*,*) 'error: singular matrix'
        END IF
    END DO

    CALL cpu_time(stop_time)
    solvetime = stop_time - start_time



    !CALL cpu_time(start_time)
    !CALL dgesv(n,m,matrix,n,ipiv,bb,n,info)
    !CALL cpu_time(stop_time)

    end subroutine solveGauss

    subroutine solveThomas(bb, matrix, n, time)
    DOUBLE PRECISION, DIMENSION(:), intent(inout)    :: bb
    DOUBLE PRECISION, DIMENSION(:,:), intent(in)  :: matrix
    INTEGER, intent(in) :: n
    DOUBLE PRECISION, intent(out) :: time

    double precision, dimension(n) :: alpha, beta
    double precision :: a, b, c
    double precision :: start_time, stop_time

    integer :: i

    call CPU_TIME(start_time)
    alpha(n) = 0d0
    beta(n) = 0d0

    do i = n, 2, -1
        a = matrix(i, i - 1)
        b = matrix(i, i)
        if (i == n) then
            c = 0d0
        else
            c = matrix(i, i +1)
        end if
        alpha(i-1) = -a / (b + c* alpha(i))
        beta(i-1) = (bb(i) - c * beta(i)) / (b + c * alpha(i))
    end do
    bb(1) = (bb(1) - matrix(1,2) * beta(1)) / (matrix(1,1) + matrix(1,2) * alpha(1))
    do i = 1, n-1
        bb(i + 1) = alpha(i) * bb(i) + beta(i)
    end do
    call CPU_TIME(stop_time)
    time = stop_time - start_time
    end subroutine solveThomas

    
    subroutine showMatrix(a, name)
    implicit none
    double precision, dimension(:,:), intent(in) :: a
    character(*), intent(in), optional :: name
    integer :: i, j
    write(*,*) "*** Show Matrix ", name, " ***"
    do i = 1, size(a,1)
        do j = 1, size(a,2)
            write(*,'(F8.4)',advance="no") a(i,j)
        end do
        write(*,*)
    end do
    end subroutine showMatrix

    END MODULE LinearSystems

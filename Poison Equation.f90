!  PoisonEquation.f90 
!
!  FUNCTIONS:
!  PoisonEquation - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: PoisonEquation
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program PoisonEquation
    use LAPACK95
    use LinearSystems
    implicit none
    double precision :: x, start_time, stop_time
    double precision :: f
    integer :: n = 700
    double precision :: delta_x 
    double precision, dimension(:,:), allocatable :: d, l , u , t, diag, subdiag
    double precision, dimension(:), allocatable :: zGauss , zLAPACK, zThomas, zThomasLAPACK
    double precision :: timelu = 0d0, timesolve = 0d0, timeGauss, timeThomas, timeGaussLapack, timeThomasLapack
    integer, dimension(:), allocatable :: ip
    integer :: i, j, k, info
    
    allocate(d(n-1,n-1))
    allocate(diag(n-1,n-1))
    allocate(subdiag(n-1,n-1))
    allocate(l(n-1,n-1))
    allocate(u(n-1,n-1))
    allocate(t(n-1,n-1))
    
    allocate(zGauss(n-1))
    allocate(zThomas(n-1))
    allocate(zThomasLAPACK(n-1))
    allocate(zLAPACK(n-1))
    
    allocate(ip(n-1))
    d = 0d0
    l = 0d0
    u = 0d0
    t = 0d0
    diag = 0d0
    subdiag = 0d0
    
    zGauss = 0d0
    zThomas = 0d0
    zThomasLAPACK = 0d0
    zLAPACK = 0d0
    
    delta_x = 1d0 / n
    !Bestimmung der Matrix für die Ableitung - Triagonale Matrix
    d(1,1) = 2d0
    diag(1,1) = 2d0
    d(1,2) = -1d0
    subdiag(1,2) = -1d0
    do i = 2, n - 2
        d(i, i-1) = -1d0
        d(i, i) = 2d0
        d(i, i +1) = -1d0
        diag(i,i) = 2d0
        subdiag(i, i - 1) = -1d0
        subdiag(i,i +1) = -1d0
    end do
    d(n-1,n-2) = -1d0
    d(n-1, n-1) = 2d0
    diag(n-1, n-1) =2d0
    subdiag(n-1, n-2) = -1d0
    !Bestimmung von b (im z Vektor weil er in lufak.. zu z wird)
    do i = 1, n - 1
        zGauss(i) = f(i * delta_x)
        zThomas(i) = f(i * delta_x)
        zThomasLAPACK(i) = f(i * delta_x)
        zLAPACK(i) = f(i * delta_x)
    end do
    !!Matrix zum Test anzeigen
    !do i = 1, n-1
    !    do j = 1, n-1
    !        write(*,'(F8.4)', advance="no") d(i,j)
    !    end do
    !    write(*,*) ""
    !end do
    !write(*,*)
    
    
    !GAUß
    print*, "GAUSS"
    u = d
    
    call lufaktorisierungGauss(u, l, zGauss, n-1,n-1, timelu)
    print*, "TIME (FAKTORISATION): " , timelu
            !!Matrix zum Test anzeigen
            !write(*,*) "L Matrix"
            !do i = 1, n-1
            !    do j = 1, n-1
            !        write(*,'(F8.4)', advance="no") l(i,j)
            !    end do
            !    write(*,*) ""
            !end do
            !
            !!Matrix zum Test anzeigen
            !write(*,*) "U Matrix"
            !do i = 1, n-1
            !    do j = 1, n-1
            !        write(*,'(F8.4)', advance="no") u(i,j)
            !    end do
            !    write(*,*) ""
            !end do
            !
            !write(*,*) "LU - A Matrix"
            !t = Matmul(l,u)
            !do i = 1, n-1
            !    do j = 1, n-1
            !        write(*,'(F8.4)', advance="no") t(i,j) - d(i,j)
            !    end do
            !    write(*,*) ""
            !end do
            !!Als Test die Build in Methode aufrufen (Intel Compiler)
            !write(*,*) "Intel RF"
            !call dgetrf(n-1,n-1,d,n-1,ip,k)
            !do i = 1, n-1
            !do j = 1, n-1
            !    write(*,'(F8.4)', advance="no")  d(i,j)
            !end do
            !write(*,*) ""
            !end do
    
    call solveGauss(zGauss, u, n-1,n-1,timesolve)
    print*, "TIME (SOlVE) : " , timesolve
    print*, "TIME (TOTAL) : " , timesolve + timelu
    timeGauss = timesolve + timelu
    u = d
    !LAPACK ROUTINE (vom INTEL COMPILER + GNU COMPILER)
    call CPU_TIME(start_time)
    call dgesv(n-1, 1, u, n-1,ip,zLAPACK,n-1, info)
    
    call CPU_TIME(stop_time)
    
    !LAPACK ALGORITHMUS (vom GNU COMPILER)
    !call dptsv(n-1, 1, diag, subdiag, zThomasLAPACK,n-1, info)
    print*,"LAPACK TIME (TOTAL) : " , stop_time - start_time
    timeGaussLapack = stop_time - start_time
    
    !THOMAS ALGORITHM
    print*, ""
    print*, "THOMAS"
    u = d 
    call solveThomas(zThomas, u, n-1, timeThomas)
    print*, "THOMAS TIME (TOTAL) : ", timeThomas
    
    do j = 1, n - 1
        if (abs(zGauss(j) - zLAPACK(j)) > 1d-5) then
            print*, "Fehler in Zeile bei Gauss" , j 
        end if
        if (abs(zThomas(j) - zLAPACK(j)) > 1d-5) then
            print*, "Fehler in Zeile bei Thomas" , j 
        end if
    end do
    
    
    deallocate(d)
    deallocate(l)
    deallocate(u)
    deallocate(t)
    deallocate(diag)
    deallocate(subdiag)
    
    deallocate(zGauss)
    deallocate(zThomas)
    deallocate(zThomasLAPACK)
    deallocate(zLAPACK)
    
    deallocate(ip)
    
    end program PoisonEquation

    function f(x)
        double precision :: f
        double precision, intent(in) :: x
        f = 100d0 * exp(-10d0 * x)
    end function f
    
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
    double precision :: f, exakt
    integer :: n = 700
    double precision :: delta_x 
    double precision, dimension(:,:), allocatable :: d, l , u , t
    double precision, dimension(:), allocatable :: zGauss , zLAPACK, zThomas, zThomasLAPACK, diag, subdiag
    double precision :: timelu = 0d0, timesolve = 0d0, timeGauss, timeThomas, timeGaussLapack, timeThomasLapack
    integer, dimension(:), allocatable :: ip
    integer :: i, j, k, info
    
    open(unit=100, file="time.dat",action="write")
    open(unit=101, file="lsgn10.dat",action="write")
    open(unit=102,file="lsgn100.dat", action ="write")
    
    do n = 10, 8000, 10
    
    
        allocate(d(n-1,n-1))
        allocate(l(n-1,n-1))
        allocate(u(n-1,n-1))
        allocate(t(n-1,n-1))
    
        allocate(diag(n-1))
        allocate(subdiag(n-2))
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
        
        !Bestimmung der Matrix für die Ableitung - Tridiagonale Matrix
        d(1,1) = 2d0
        diag(1) = 2d0
        d(1,2) = -1d0
        subdiag(1) = -1d0
        do i = 2, n - 2
            d(i, i-1) = -1d0
            d(i, i) = 2d0
            d(i, i +1) = -1d0
            diag(i) = 2d0
            subdiag(i) = -1d0
        end do
        d(n-1,n-2) = -1d0
        d(n-1, n-1) = 2d0
        diag(n-1) =2d0
        
        
        
        !Bestimmung von b (im z Vektor weil er in lufak.. zu z wird)
        do i = 1, n - 1
            zGauss(i) = delta_x **2 * f(i * delta_x)
            zThomas(i) = delta_x **2 *f(i * delta_x)
            zThomasLAPACK(i) = delta_x **2 *f(i * delta_x)
            zLAPACK(i) = delta_x **2 *f(i * delta_x)
        end do
    
    
        !GAUß
        if (n < 600) then
            u = d
            call lufaktorisierungGauss(u, l, zGauss, n-1,n-1, timelu)
            !print*, "TIME (FAKTORISATION): " , timelu
                    !!Matrix zum Test anzeigen
                    !call showMatrix(l,"L")
                    !call showMatrix(u,"u")
                    !t = Matmul(l,u)
                    !call showMatrix(t,"LU - A Matrix")
                    !!Als Test die Build in Methode aufrufen (Intel Compiler)
                    !call dgetrf(n-1,n-1,d,n-1,ip,k)
                    !call showMatrix(d,"INTEL RL Test")
    
            call solveGauss(zGauss, u, n-1,n-1,timesolve)
            timeGauss = timesolve + timelu
        end if
        
        
        u = d
        !LAPACK ROUTINE (vom INTEL COMPILER + GNU COMPILER)
        call CPU_TIME(start_time)
        call dgesv(n-1, 1, u, n-1,ip,zLAPACK,n-1, info)  
        call CPU_TIME(stop_time)
        timeGaussLapack = stop_time - start_time
    !LAPACK ALGORITHMUS (vom GNU COMPILER)
    !call CPU_TIME(start_time)
    !call dptsv(n-1, 1, diag, subdiag, zThomasLAPACK,n-1, info)
    !call CPU_TIME(stop_time)
    !timeThomasLAPACK = stop_time - start_time
        
    
        !THOMAS ALGORITHMUS
        u = d 
        call solveThomas(zThomas, u, n-1, timeThomas)
        
        
        !Vergleiche die Ergebnisse
        do j = 1, n - 1
            if (n < 600) then
                if (abs(zGauss(j) - zLAPACK(j)) > 1d-5) then
                    print*, "Fehler in Zeile bei Gauss" , j 
                end if
            end if
            if (abs(zThomas(j) - zLAPACK(j)) > 1d-5) then
                print*, "Fehler in Zeile bei Thomas Gauss" , j 
            end if
!if (abs(zThomas(j) - zThomasLAPACK(j)) > 1d-5) then
!    print*, "Fehler in Zeile bei Thomas" , j 
!end if
        end do
    
        !Ausgabe für n=10 und n= 100 in die Datein schreiben
        if (n == 10) then
            write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), 0d0, 0d0, 0d0, 0d0, 0d0 !Randbed.
!write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8.4)'), 0d0, 0d0, 0d0, 0d0, 0d0, 0d0!Randbed.
            do i = 1, n - 1
                x = i * delta_x
                write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), x, exakt(x), zGauss(i), zThomas(i), zLAPACK(i)
!write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8.4)'), x, exakt(x), zGauss(i), zThomas(i), zLAPACK(i), zThomasLAPACK
            end do       
            write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), 1d0, 0d0, 0d0, 0d0, 0d0
!write(101, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8,4)'), 1d0, 0d0, 0d0, 0d0, 0d0, 0d0
        else if (n == 100) then
             write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), 0d0, 0d0, 0d0, 0d0, 0d0 !Randbed.
!write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8.4)'), 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 !Randbed.
            do i = 1, n - 1
                x = i * delta_x
                write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), x, exakt(x), zGauss(i), zThomas(i), zLAPACK(i)
!write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8.4)'), x, exakt(x), zGauss(i), zThomas(i), zLAPACK(i), zThomasLAPACK(i)
            end do       
            write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4)'), 1d0, 0d0, 0d0, 0d0, 0d0
!write(102, '(F8.4 F8.4 F8.4 F8.4 F8.4 F8.4)'), 1d0, 0d0, 0d0, 0d0, 0d0, 0d0
        end if
        
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
        if (n < 600) then
            write(100, '(I8.3 F16.6 F16.6 F16.6)') n,  timeGaussLAPACK, timeThomas, timeGauss
!write(100, '(I8.4 F8.4 F8.4 F8.4 F8.4') n, timeGaussLAPACK, timeThomas, timeThomasLAPACK, timeGauss, 
        else     
            write(100, '(I8.3 F16.6 F16.6 F16.6)') n,  timeGaussLAPACK, timeThomas
!write(100, '(I8.4 F8.4 F8.4 F8.4 F8.4') n, timeGaussLAPACK, timeThomas, timeThomasLAPACK
        end if
        
    end do
    
    close(100)
    close(101)
    close(102)
    end program PoisonEquation

    function f(x)
        double precision :: f
        double precision, intent(in) :: x
        f =  100d0 * exp(-10d0 * x)
    end function f
    function exakt(x)
        double precision :: exakt
        double precision, intent(in) :: x
        exakt = 1d0 - (1d0 - exp(-10d0)) * x - exp(-10d0 * x)
    end function
    
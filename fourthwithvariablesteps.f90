program oneToTwelveNavier
!This is a fortran implementation of Burger's Equation
!FORTRAN implementaion of 12 Steps to Navier Stokes by Lorena Barba(*****https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/*****)  
implicit none
integer::n
integer::i,tsim,imv,j,k
real,allocatable, dimension(:)::x,xtpo
real::dx,dt,l,cou,t,xa,ta,pi,nu,lo
pi=4*atan(1.000000)
write(*,*)"Please enter number of divisions in the domain,(atleast 20)"
read *,n
write(*,*)"Please enter total length"
read *,l
write(*,*)"Please enter simulation time in s"
read *,t
write(*,*)"Please enter viscosity nu in m/s"
read *,nu
allocate(x(n),xtpo(n))
dx=l/n
dt=dx*nu
!cou=u*(dt/dx)
tsim=int(t/dt)
k=int(n/2)
!	if (cou>=1) then
!		write (*,*)"This will blow up, why?"
!		write (*,*)"courant number is"
!		write (*,*) cou
!	else
!		write (*,*)"continue"
!		write (*,*)"courant number is"
!		write (*,*) cou
!	end if
open(1,file="InitialFunc.dat", status='replace')
open(2,file="NumericalSolution.dat", status='replace' )
open(3,file="result.plt", status='replace')

!initial file sawtooth
	do i=0,n
        if ((i<=k)) then
            x(i)=4+((3*(i*dx))/(k*dx))      
        else
            x(i)=1+((3*((i-k)*dx))/((n-k)*dx)) 
        end if
		write(1,*)i*dx,x(i)
	end do
!processing
    do j=0,tsim
		do i=1,n
            x(1)=x(n)
			xtpo(i)=x(i)-(((dt*x(i))/dx)*(x(i)-x(i-1)))+(((nu*dt)/(dx**2))*(x(i+1)-(2*x(i))+x(i-1)))
            write(*,*) x(i)		
        end do
		x=xtpo
	end do
do i=1,n
write(2,*)i*dx,x(i)
end do
write(3,*)'set xlabel "spatial dimension"'
write(3,*)'set ylabel "Value"'
write(3,*)'set title "Burgers Equation"'
write(3,*)'plot "InitialFunc.dat" with line, "NumericalSolution.dat" with line'
CALL SYSTEM('gnuplot -p result.plt')
close(1)
close(2)
close(3)
end program oneToTwelveNavier

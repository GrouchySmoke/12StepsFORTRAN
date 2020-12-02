program oneToTwelveNavier
!This is a fortran implementation of 1D non linear convection
!FORTRAN implementaion of 12 Steps to Navier Stokes by Lorena Barba(*****https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/*****)  
implicit none
integer::n,a
integer::i,tsim,j
real,dimension(:),allocatable::T,Ttpo,x
real::dx,dt,l,cou,t1,u,alpha
real,parameter::pi=4*atan(1.0d0)
write(*,*)"Please enter number of divisions in the domain,(atleast 20)"
read *,n
write(*,*)"Please enter total length"
read *,l
write(*,*)"Please enter time step value in s"
read *,dt
write(*,*)"Please enter simulation time in s"
read *,t1
write(*,*)"Please enter thermal conductivity"
read *,alpha
allocate(T(n),Ttpo(n),x(n))
dx=l/n
u=2
cou=alpha*(dt/dx**2)
tsim=int(t1/dt)
if (cou>=0.500000) then
write (*,*)"This will blow up, why?"
else
write (*,*)"continue"
end if
!initialising a step function
	do i=0,n
        x(i)=i*dx!spatial array
        T(i)=500*sin(pi*(x(i)/l))!use different functions and modify
	end do
open(1,file="attzero.dat", status='replace')
open(2,file="attfinal.dat", status='replace')
open(3,file="attzero.plt", status='replace')
open(4,file="attfinal.plt", status='replace')
!writing to file at initial time step
	do i=0,n
		write(1,*)(i*dx),T(i)
	end do
!processing the 1D linear convection equation
	do j=1,tsim
		do i=1,n-1
			Ttpo(i)=T(i)+(cou*(T(i+1)-(2*T(i))+T(i-1)))
            !write (*,*) Ttpo(i)
		end do
        !write (*,*)
		T=Ttpo
	end do
!writing final values	
	do i=0,n
		write(2,*)(i*dx),T(i)
	end do
write(3,*)'set xlabel "spatial dimension"'
write(3,*)'set ylabel "Value"'
write(3,*)'set title "At time=0"'
write(3,*)'plot "attzero.dat" with line,"attfinal.dat" with line'
CALL SYSTEM('gnuplot -p attzero.plt')
close(1)
close(2)
close(3)
close(4)
!deallocate(x,xtpo)
!Deallocating throws a SIGABRT signal when compiled in gfortran
end program oneToTwelveNavier

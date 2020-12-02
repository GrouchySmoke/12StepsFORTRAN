program oneToTwelveNavier
!This is a fortran implementation of 1D non linear convection
!FORTRAN implementaion of 12 Steps to Navier Stokes by Lorena Barba(*****https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/*****)  
implicit none
integer::n,a
integer::i,tsim,imv,j
real,dimension(:),allocatable::x,xtpo
real::dx,dt,l,cou,t,u
write(*,*)"Please enter number of divisions in the domain,(atleast 20)"
read *,n
imv=n+1
write(*,*)"Please enter total length"
read *,l
write(*,*)"Please enter time step value in s"
read *,dt
write(*,*)"Please enter simulation time in s"
read *,t
!write(*,*)"Please enter wave speed c in m/s"
!read *,c
allocate(x(imv),xtpo(imv))
dx=l/n
u=2
cou=u*(dt/dx)
tsim=int(t/dt)
if (cou>=1) then
write (*,*)"This will blow up, why?"
write (*,*)"courant number is"
write (*,*) cou
else
write (*,*)"continue"
write (*,*)"courant number is"
write (*,*) cou
end if
!initialising a step function
	do i=0,n
		if ((i>5).and.(i<10)) then
			x(i)=u
		else 
			x(i)=1
		end if
	end do
open(1,file="attzero.dat", status='replace')
open(2,file="attfinal.dat", status='replace')
open(3,file="attzero.plt", status='replace')
open(4,file="attfinal.plt", status='replace')
!writing to file at initial time step
	do i=0,n
		write(1,*)(i*dx),x(i)
	end do
!processing the 1D linear convection equation
	do j=0,tsim
		do i=1,n
			xtpo(i)=x(i)-(x(i)*(dt/dx))*(x(i)-x(i-1))
		end do
		x=xtpo
	end do
!writing final values	
	do i=0,n
		write(2,*)(i*dx),x(i)
	end do
write(3,*)'set xlabel "spatial dimension"'
write(3,*)'set ylabel "Value"'
write(3,*)'set title "At time=0"'
write(3,*)'plot "attzero.dat" with line'
CALL SYSTEM('gnuplot -p attzero.plt')
write(4,*)'set xlabel "spatial dimension"'
write(4,*)'set ylabel "Value"'
write(4,*)'set title "At time=final"'
write(4,*)'plot "attfinal.dat" with line'
CALL SYSTEM('gnuplot -p attfinal.plt')
close(1)
close(2)
close(3)
close(4)
!read *,a
!deallocate(x,xtpo)
!Deallocating throws a SIGABRT signal when compiled in gfortran
end program oneToTwelveNavier

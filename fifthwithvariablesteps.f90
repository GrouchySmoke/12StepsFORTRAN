program oneToTwelveNavier
!This is a fortran implementation of 1D linear advection 
!FORTRAN implementaion of 12 Steps to Navier Stokes by Lorena Barba(*****https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/*****)  
implicit none
integer::nx,ny,tt
integer::i,tsim,imvx,imvy,j
real,allocatable, dimension(:,:)::x,xtpo
real::dx,dy,dt,lx,ly,cou,cou1,cou2,c,t
write(*,*)"Please enter number of divisions in the X domain,(atleast 20)"
read *,nx
write(*,*)"Please enter number of divisions in the Y domain,(atleast 20)"
read *,ny
imvx=nx-1
imvy=ny-1
write(*,*)"Please enter total X length"
read *,lx
write(*,*)"Please enter total y length"
read *,ly
write(*,*)"Please enter time step value in s"
read *,dt
write(*,*)"Please enter simulation time in s"
read *,t
write(*,*)"Please enter wave speed c in m/s"
read *,c
allocate(x(nx,ny),xtpo(nx,ny))
dx=lx/nx
dy=ly/ny
cou1=c*(dt/dx)
cou2=c*(dt/dy)
cou=max(cou1,cou2)
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
	do i=0,nx
        do j=0,ny
		    if (((i>5).and.(i<10)).and.((j>5).and.(j<10))) then
			    x(i,j)=5
		    else 
			    x(i,j)=0
		    end if
        end do
	end do
open(1,file="attzero.dat", status='replace')
open(2,file="attfinal.dat", status='replace')
open(3,file="attzero.plt", status='replace')
open(4,file="attfinal.plt", status='replace')
!writing to file at initial time step
	do i=0,nx
        do j=0,ny
		    write(1,*) i,j,x(i,j)
        end do
	end do
!processing the 1D linear convection equation
	do tt=0,tsim
		do i=1,nx
            do j=1,ny
			    xtpo(i,j)=x(i,j)-cou1*(x(i,j)-x(i-1,j))-cou2*(x(i,j)-x(i,j-1))
            end do
		end do
		x=xtpo
	end do
!writing final values	
	do i=0,nx
        do j=0,ny
		    write(2,*)i,j,x(i,j)
        end do
	end do
close(1)
close(2)
write(3,*)'set xlabel "X dimension"'
write(3,*)'set ylabel "Y dimension"'
write(3,*)'set zlabel "Value"'
write(3,*)'set title "At time=0"'
write(3,*)'set hidden3d'
write(3,*)'splot "attzero.dat"'
CALL SYSTEM('gnuplot -p attzero.plt')
write(4,*)'set xlabel "X dimension"'
write(4,*)'set ylabel "Y dimesion"'
write(4,*)'set ylabel "Value"'
write(4,*)'set title "At time=final"'
write(3,*)'set hidden3d'
write(4,*)'splot "attfinal.dat"'
CALL SYSTEM('gnuplot -p attfinal.plt')
close(3)
close(4)
deallocate(x,xtpo)
end program oneToTwelveNavier
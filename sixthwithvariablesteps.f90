program oneToTwelveNavier
!This is a fortran implementation of 2D linear convection 
!FORTRAN implementaion of 12 Steps to Navier Stokes by Lorena Barba(*****https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/*****)  
implicit none
integer::nx,ny,tt
integer::i,tsim,imvx,imvy,j
real,allocatable, dimension(:,:)::u,utpo,v,vtpo
real::dx,dy,dt,lx,ly,cou,cou1,cou2,c,t
real::start,end
write(*,*)"Please enter number of divisions in the X domain,(atleast 20)"
read *,nx
write(*,*)"Please enter number of divisions in the Y domain,(atleast 20)"
read *,ny
write(*,*)"Please enter total X length"
read *,lx
write(*,*)"Please enter total y length"
read *,ly
write(*,*)"Please enter time step value in s"
read *,dt
write(*,*)"Please enter simulation time in s"
read *,t
!write(*,*)"Please enter wave speed c in m/s"
!read *,c
allocate(u(nx,ny),utpo(nx,ny),v(nx,ny),vtpo(nx,ny))
dx=lx/nx
dy=ly/ny
write(*,*)"allocate worked"
!cou1=u*(dt/dx)
!cou2=v*(dt/dy)
!cou=max(cou1,cou2)
tsim=int(t/dt)
!	if (cou>=1) then
	!	write (*,*)"This will blow up, why?"
		!write (*,*)"courant number is"
		!write (*,*) cou
	!else
		!write (*,*)"continue"
		!write (*,*)"courant number is"
		!write (*,*) cou
	!end if
!initialising a step function in u and v
	do i=0,nx
        do j=0,ny
		    if (((i>5).and.(i<10)).and.((j>5).and.(j<10))) then
			    u(i,j)=5
		    else 
			    u(i,j)=1
		    end if
        end do
	end do
	do i=0,nx
        do j=0,ny
		    if (((i>5).and.(i<10)).and.((j>5).and.(j<10))) then
			    v(i,j)=5
		    else 
			    v(i,j)=1
		    end if
        end do
	end do
write(*,*)"initialisation worked"
open(1,file="attzero.dat", status='replace')
open(2,file="attfinalu.dat", status='replace')
open(3,file="attfinalv.dat", status='replace')
open(4,file="attzero.plt", status='replace')
open(5,file="attfinalu.plt", status='replace')
open(6,file="attfinalv.plt", status='replace')
!writing to file at initial time step
	do i=0,nx
        do j=0,ny
		    write(1,*) i,j,sqrt((u(i,j)**2)+(v(i,j)**2))
        end do
	end do
!processing the 2D linear convection equation
call cpu_time(start)
	do tt=0,tsim
		do i=1,nx
            do j=1,ny
			    utpo(i,j)=u(i,j)-(u(i,j)*(dt/dx)*(u(i,j)-u(i-1,j)))-(v(i,j)*(dt/dx)*(u(i,j)-u(i,j-1)))
                vtpo(i,j)=v(i,j)-(u(i,j)*(dt/dx)*(v(i,j)-v(i-1,j)))-(v(i,j)*(dt/dx)*(v(i,j)-v(i,j-1)))
            end do
		end do
		u=utpo
        v=vtpo
	end do
call cpu_time(end)
!write(*,*)"time of execution is", end-start, "s"
!print *,"time of execution is (f20.15) s",end-start
!writing final values	
	do i=0,nx
        do j=0,ny
		    write(2,*)i,j,u(i,j)
            write(3,*)i,j,v(i,j)
        end do
	end do

write(4,*)'set xlabel "X dimension"'
write(4,*)'set ylabel "Y dimension"'
write(4,*)'set zlabel "Value of velocity"'
write(4,*)'set title "At time=0"'
write(4,*)'set hidden3d'
write(4,*)'splot "attzero.dat"'
CALL SYSTEM('gnuplot -p attzero.plt')
write(5,*)'set xlabel "X dimension"'
write(5,*)'set ylabel "Y dimesion"'
write(5,*)'set zlabel "Value of u velocity"'
write(5,*)'set title "At time=final"'
write(5,*)'set hidden3d'
write(5,*)'splot "attfinalu.dat"'
CALL SYSTEM('gnuplot -p attfinalu.plt')
write(6,*)'set xlabel "X dimension"'
write(6,*)'set ylabel "Y dimesion"'
write(6,*)'set zlabel "Value of v velocity"'
write(6,*)'set title "At time=final"'
write(6,*)'set hidden3d'
write(6,*)'splot "attfinalv.dat"'
CALL SYSTEM('gnuplot -p attfinalv.plt')
close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
deallocate(u,utpo)
deallocate(v,vtpo)
end program oneToTwelveNavier
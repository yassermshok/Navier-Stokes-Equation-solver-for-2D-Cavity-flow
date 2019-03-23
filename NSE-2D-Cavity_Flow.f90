module Vars
    
implicit none

!! i, j , and k are index variables for x, y and z
!! it is index for time steps
!! nx, ny adre number of cells in x and why directions
!! dx and dy are the cell size in x and y directions
!! dt is the time step
!! Lx and Ly are the length of the domain in x and y directions
!! P is an array to store Pressure values
!! PP is an array to store pressure values of previous iteration
!! B is an array for source terms
!! frame is an index for different plot frames
!! iter is index for iteration number


integer(8) i, j, k, nx, ny, frame, iter, plot_frame, it
real(8) dx, dy, Lx, Ly, dt, rho, neu
real(8), allocatable :: P(:,:), B(:,:), PP(:,:), u(:,:), v(:,:), un(:,:), vn(:,:)

end module Vars

program NSE
    
    use Vars
    
    implicit none
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Initialization**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!! Domain Dimensions
    
    Lx = 2
    Ly = 2
    
!! Number of points

    nx = 20
    ny = 20
    
!! Computing dx, dy, and setting dt

    dx = Lx/nx
    dy = Ly/ny
    dt = 0.01

!! Fluid Properties
    
    rho = 1.2
    neu = 0.1
    
!! Intializing flow variables
    
    allocate ( P(0:nx, 0:ny), PP(0:nx, 0:ny), B(0:nx, 0:ny), u(0:nx, 0:ny), v(0:nx, 0:ny), un(0:nx, 0:ny), vn(0:nx, 0:ny) )

    P = 0
    PP = 0
    B = 0
    un = 0
    vn = 0
    un(0:nx, ny) = 1
    
!! Writting a plot file with the intial conditions

    plot_frame = 0
    frame = 1000
    call write_vtk
    print*, "Completed Initialization" 
        
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Solution**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!! outer loop for time steps.

do it = 1, 500
    
    u = un
    v = vn
    
!! Nested Loop to calculate B(i,j) in Poisson equation

    do i = 1, nx-1
    do j  = 1, ny-1
      B(i,j) = -rho * ( ((u(i+1,j)-u(i-1,j))/2/dx)**2 + ((v(i,j+1)-v(i,j-1))/2/dy)**2 +2*(u(i,j+1)-u(i,j-1))/2/dy*(v(i+1,j)-v(i-1,j))/2/dx ) + rho/dt*( ((u(i+1,j)-u(i-1,j))/2/dx) + ((v(i,j+1)-v(i,j-1))/2/dy) )
    end do
    end do
 
!! Iterating to calculate Pressure

do iter = 1, 500

PP = P 
     do i = 1, nx-1
     do j = 1, ny-1
       P(i,j) = (dy**2 *( PP(i-1,j) + PP(i+1,j) ) + dx**2 *( PP(i,j-1) + PP(i,j+1) )-B(i,j)*dx**2*dy**2 ) /(2*(dx**2 + dy**2))
     end do
     end do
     
!! Imposing Pressure Boundary conditions

     do i = 0, nx
         P(i, ny) = 0
         P(i,0) =  P(i,1)
     end do 

      do j = 0, ny
         P(nx,j) = P(nx-1,j)
         P(0,j) =  P(1,j)
     end do 
end do


!! Calculating new time step velcoity

do i = 1, nx-1
do j = 1, ny-1
        
un(i,j) = dt *( ( -1/rho*(P(i+1,j)-P(i-1,j))/2/dx + neu*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx/dx+(u(i,j+1)-2*u(i,j)+u(i,j-1))/dy/dy)) - u(i,j)*(u(i,j)-u(i-1,j))/dx  - v(i,j)*(u(i,j)-u(i,j-1))/dy   ) + u(i,j)

vn(i,j) = dt *( ( -1/rho*(P(i,j+1)-P(i,j-1))/2/dy + neu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx/dx+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy/dy)) - u(i,j)*(v(i,j)-v(i-1,j))/dx  - v(i,j)*(v(i,j)-v(i,j-1))/dy   ) + v(i,j)

end do
end do

!! writing a plot file for the current time step.

plot_frame = plot_frame+1
if (plot_frame .eq. 1 ) then
plot_frame = 0
print*, "Completed time step no. ",  it  
frame = frame +1
call write_vtk
end if

end do

end program NSE



subroutine write_vtk
    use Vars
    implicit none
    
    character (len = 1024) ::  plot

    write(plot,"(A7,I4,A4)")  "ABYtest", frame, ".vtk"
    open(unit=31, file=plot)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing file header and geometry**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
   
   write (31, 501)
   write (31, 502) 
   write (31, 503) 
   write (31, 504) 
   write (31, 505) (nx+1)*(ny+1)*2
   
   do k = 0,1
     do i = 0, nx
        do j = 0, ny
             
          write (31,601) i*dx, j*dy, float(k)
          
        end do
     end do 
   end do

    write (31,506) nx*ny,  nx*ny*9

    do i = 1, nx
    do j = 1, ny
    write (31,602) ((ny+1)*(i-1) + j)-1, ((ny+1)*(i)+j)-1, ((ny+1)*(i) + j+1)-1, ((ny+1)*(i-1) + j+1)-1, (nx+1)*(ny+1)+((ny+1)*(i-1) + j)-1, (nx+1)*(ny+1)+((ny+1)*(i) + j)-1,(nx+1)*(ny+1)+((ny+1)*(i) + j+1)-1,(nx+1)*(ny+1)+((ny+1)*(i-1) + j+1)-1        
    end do
    end do
    
    write (31, 507) nx*ny
    
    do i = 1, nx*ny
    write (31,603) 12
    end do
    
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing Values of Pressure**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
    
    write (31, 508)  (nx+1)*(ny+1)*2
    write (31, 509)
    write (31, 510)
    
    do k = 0, 1  
    do i = 0, nx
    do j = 0, ny
      write (31,604) P(i,j)
    end do
    end do 
    end do
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing Values of Velocity**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

    write (31, 511)
    write (31, 510)
 
    do k = 0, 1  
    do i = 0, nx
    do j = 0, ny
      write (31,604) un(i,j)
    end do
    end do 
    end do
      
      
    write (31, 512)
    write (31, 510)
    
     do k = 0, 1  
     do i = 0, nx
     do j = 0, ny
       write (31,604) vn(i,j)
    end do
    end do 
    end do
    
  write (31, 514)
  
    do k = 0, 1  
    do i = 0, nx
    do j = 0, ny
       write (31,'(3f12.8)') un(i,j), vn(i,j), 0
    end do
    end do 
    end do
   

   
501 format("# vtk DataFile Version 2.0" )   
502 format("Unstructured Grid Example")   
503 format("ASCII")       
504 format("DATASET UNSTRUCTURED_GRID")       
505 format("POINTS ", i10, " float")   
506 format("CELLS ", i10, " ", i10)  
507 format("CELL_TYPES ", i10)
508 format("POINT_DATA ", i7)
509 format("SCALARS Pressure_(Pa) float 1")
510 format("LOOKUP_TABLE default")
511 format("SCALARS u_(m/sec) float 1")
512 format("SCALARS v_(m/sec) float 1")
513 format("SCALARS w_(m/sec) float 1")    
514 format("VECTORS Velocity_(m/sec)  float")
515 format("CELL_DATA ", i7)
516 format("SCALARS Pressure_(Cell) float 1")
517 format("SCALARS u_(Cell) float 1")
518 format("SCALARS v_(Cell) float 1")
519 format("SCALARS w_(Cell) float 1")    
520 format("VECTORS Velocity_(Cell)  float")  
    
521 format("SCALARS v_(star) float 1")
522 format("SCALARS v_(corr) float 1")
523 format("SCALARS p_(corr) float 1")    
    
    
601 format(f10.5, " ",f10.5, " ",f10.5)
602 format("8 ", i7, " ",i7, " ",i7, " ",i7," ",i7," ",i7," ",i7," ",i7)
603 format(i2)
604 format(f20.10)
 
 end subroutine write_vtk

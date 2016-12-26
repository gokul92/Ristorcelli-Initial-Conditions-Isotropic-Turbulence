Program P_to_S
  implicit none                                                     !sizex,sizey,sizez are the dimension of the total box generated using code9
  integer,parameter :: sizex=96,sizey=96,sizez=96,yprocs=2,zprocs=2,nprocs=yprocs*zprocs 
  integer i,j,k,l,ierr,num1,num2,num3,num4
  real*8 fz,lz,fy,ly
  double precision,dimension(1:sizex,1:sizey,1:sizez) :: x_un,y_un,z_un,u,v,w,den,p,T
  character(len=40) :: filename
  
  !allocate(x_un(sizex,sizey,sizez),y_un(sizex,sizey,sizez),z_un(sizex,sizey,sizez),u(sizex,sizey,sizez),v(sizex,sizey,sizez),w(sizex,sizey,sizez),den(sizex,sizey,sizez),p(sizex,sizey,sizez))
  open(unit=50,file='serial_initial2.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
  !read(50,*)
  !read(50,*)
  do k=1,sizez
     do j=1,sizey
        do i=1,sizex      
           read(50,'(6E30.18)') u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),den(i,j,k),T(i,j,k)
        end do
     end do
  end do
  close(50)
  do l=0,nprocs-1
     num1=l/10
     num2=mod(l,10)
     num3=num2/10
     num4=mod(num2,10)
!!$     num3=l/10
!!$     num4=mod(l,10)
     filename='Data/'//'p'//char(num1+48)//char(num2+48)//'.dat'
!!$     filename='Turb_inflow'//char(num3+48)//char(num4+48)//'.dat'
     fz=mod(l,zprocs)*sizez/zprocs+1
     lz=fz+sizez/zprocs-1
     fy=(l/zprocs)*sizey/yprocs+1
     ly=fy+sizey/yprocs-1

     open(unit=50,file=filename,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
     write(50,*) 'variables="u","v","w","p","ro","T"'
     write(50,*) 'zone I=',sizex,'J=',sizey/yprocs,'K=',sizez/zprocs,' F=POINT'
     do k=fz,lz
        do j=fy,ly
           do i=1,sizex      
              write(50,'(6E30.18)') u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),den(i,j,k),T(i,j,k)
              !read(50,'(1x,6E30.18)') x_un(i,j,k),y_un(i,j,k),z_un(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k) 
           end do
        end do
     end do
     close(50)
  end do
  !deallocate(x_un,y_un,z_un,u,v,w,den,p)
End program P_to_S
   

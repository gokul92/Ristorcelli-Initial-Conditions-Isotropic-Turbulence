module ns

  implicit none

  include "fftw3.f"

contains

  subroutine wrap_around(f, n)

    integer n, i
    double complex, dimension(0:n-1) :: f
    double complex :: temp

    do i = 0, n/2 - 1
       temp = f(i)
       f(i) = f(n/2+i)
       f(n/2+i) = temp
    end do

  end subroutine wrap_around

  subroutine wrap_around_3d(f, nx, ny, nz)

    integer :: nx, ny, nz, i, j, k
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: f
    double complex, dimension(0:nx-1) :: input_x
    double complex, dimension(0:ny-1) :: input_y
    double complex, dimension(0:nz-1) :: input_z
    double complex :: temp

    do k = 0, nz-1
       do j = 0, ny-1
          input_x = f(:,j,k)
          call wrap_around(input_x, nx)
          f(:,j,k) = input_x
       end do
    end do

    do k = 0, nz-1
       do i = 0, nx-1
          input_y = f(i,:,k)
          call wrap_around(input_y, ny)
          f(i,:,k) = input_y
       end do
    end do

    do j = 0, ny-1
       do i = 0, nx-1
          input_z = f(i,j,:)
          call wrap_around(input_z, nz)
          f(i,j,:) = input_z
       end do
    end do

  end subroutine wrap_around_3d

  subroutine fourier_derivative(n, array)						! subroutine calculates derivative
    ! in fourier space. 
    implicit none

    integer :: n, i
    double precision :: n1
    double complex, dimension(0:n-1) :: array

    do i = 0, n-1	 							! calculating derivative of input				
       if (i .le. n/2 - 1) then 					! function. This is done bymu(i, j, k)ltiplying
          n1 = 1.0d0*i
          array(i) = array(i)*(0.0d0,1.0d0)*n1			! i*k (i is the imaginary unit) to every
       else								! element in the array. Here, k runs from
          n1 = 1.0d0*(i-n)
          array(i) = array(i)*(0.0d0,1.0d0)*n1			! -l/2 to l/2 - 1. Note that the
       end if								! array is in wrapped around format only.
    end do

    array(n/2) = 0.0d0							! setting oddball wavenumber = 0

  end subroutine fourier_derivative

  subroutine x_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)

    implicit none

    integer :: i, j, k, nx, ny, nz		
    double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
    double complex, dimension(0:nx-1) :: input_x
    integer*8 ::plan_x_fwd, plan_x_bwd
    integer :: fwd, bwd, deriv

    call dfftw_plan_dft_1d(plan_x_fwd, nx, input_x, input_x, fftw_forward, fftw_measure)
    call dfftw_plan_dft_1d(plan_x_bwd, nx, input_x, input_x, fftw_backward, fftw_measure)

    if (fwd .eq. 1) then
       do k = 0, nz-1
          do j = 0, ny-1										! transforming to fourier space.
             input_x = input_3d(:, j, k)							! copy to 1d array from 3d array.
             call dfftw_execute_dft(plan_x_fwd, input_x, input_x)				! call fftw 1d fft subroutine.
             if (deriv .eq. 1) then
                call fourier_derivative(nx, input_x)					! call derivative subroutine.
             end if
             input_3d(:, j, k) = input_x							! copy back to 3d input array from
          end do											! 1d array.
       end do
       input_3d = input_3d/nx										! renormalization after x-fft.
    end if

    if (bwd .eq. 1) then
       do k = 0, nz-1
          do j = 0, ny-1										! transforming back to physical space.
             input_x = input_3d(:, j, k)
             call dfftw_execute_dft(plan_x_bwd, input_x, input_x)
             input_3d(:, j, k) = input_x
          end do
       end do
    end if

    call dfftw_destroy_plan(plan_x_fwd)
    call dfftw_destroy_plan(plan_x_bwd)

  end subroutine x_derivative

  subroutine y_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)

    implicit none

    integer :: i, j, k, nx, ny, nz		
    double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
    double complex, dimension(0:ny-1) :: input_y
    integer*8 ::plan_y_fwd, plan_y_bwd
    integer :: fwd, bwd, deriv

    call dfftw_plan_dft_1d(plan_y_fwd, ny, input_y, input_y, fftw_forward, fftw_measure)
    call dfftw_plan_dft_1d(plan_y_bwd, ny, input_y, input_y, fftw_backward, fftw_measure)

    if (fwd .eq. 1) then
       do k = 0, nz-1
          do i = 0, nx-1
             input_y = input_3d(i, :, k)
             call dfftw_execute_dft(plan_y_fwd, input_y, input_y)
             if (deriv .eq. 1) then
                call fourier_derivative(ny, input_y)
             end if
             input_3d(i, :, k) = input_y
          end do
       end do
       input_3d = input_3d/ny
    end if

    if (bwd .eq. 1) then
       do k = 0, nz-1
          do i = 0, nx-1
             input_y = input_3d(i, :, k)
             call dfftw_execute_dft(plan_y_bwd, input_y, input_y)
             input_3d(i, :, k) = input_y
          end do
       end do
    end if

    call dfftw_destroy_plan(plan_y_fwd)
    call dfftw_destroy_plan(plan_y_bwd)

  end subroutine y_derivative

  subroutine z_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)

    implicit none

    integer :: i, j, k, nx, ny, nz		
    double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
    double complex, dimension(0:nz-1) :: input_z
    integer*8 ::plan_z_fwd, plan_z_bwd
    integer :: fwd, bwd, deriv

    call dfftw_plan_dft_1d(plan_z_fwd, nz, input_z, input_z, fftw_forward, fftw_measure)
    call dfftw_plan_dft_1d(plan_z_bwd, nz, input_z, input_z, fftw_backward, fftw_measure)

    if (fwd .eq. 1) then 
       do j = 0, ny-1
          do i = 0, nx-1
             input_z = input_3d(i, j, :)
             call dfftw_execute_dft(plan_z_fwd, input_z, input_z)
             if (deriv .eq. 1) then
                call fourier_derivative(nz, input_z)
             end if
             input_3d(i, j, :) = input_z
          end do
       end do
       input_3d = input_3d/nz
    end if

    if (bwd .eq. 1) then
       do j = 0, ny-1
          do i = 0, nx-1
             input_z = input_3d(i, j, :)
             call dfftw_execute_dft(plan_z_bwd, input_z, input_z)
             input_3d(i, j, :) = input_z
          end do
       end do
    end if

    call dfftw_destroy_plan(plan_z_fwd)
    call dfftw_destroy_plan(plan_z_bwd)

  end subroutine z_derivative

  subroutine vdotdel(nx, ny, nz, v1, v2, v3, p1, total, deriv)

    integer :: nx, ny, nz
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: v1, v2, v3, temp, total, p1
    integer :: deriv 

    total = 0.0

    if (deriv .eq. 1) then
       temp = v1
       call x_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v1
       total=total+temp

       temp = v1
       call y_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v2
       total=total+temp

       temp = v1
       call z_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v3
       total=total+temp
    else if (deriv .eq. 2) then
       temp = v2
       call x_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v1
       total=total+temp

       temp = v2
       call y_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v2
       total=total+temp

       temp = v2
       call z_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v3
       total=total+temp
    else if (deriv .eq. 3) then
       temp = v3
       call x_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v1
       total=total+temp

       temp = v3
       call y_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v2
       total=total+temp

       temp = v3
       call z_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v3
       total=total+temp
    else if (deriv .eq. 4) then
       temp = p1
       call x_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v1
       total=total+temp

       temp = p1
       call y_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v2
       total=total+temp

       temp = p1
       call z_derivative(temp, nx, ny, nz, 1, 1, 1)
       temp = temp*v3
       total=total+temp
    end if
  end subroutine vdotdel

  subroutine poisson(input, output, nx, ny, nz)

    integer :: nx, ny, nz, i, j, k
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: input, output

    call x_derivative(input, nx, ny, nz, 1, 0, 0)
    call y_derivative(input, nx, ny, nz, 1, 0, 0)
    call z_derivative(input, nx, ny, nz, 1, 0, 0)
    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx-1
             if ( i .eq. 0 .and. j .eq. 0 .and. k .eq. 0) then
                output(i,j,k) = 0
             else if (i .le. nx/2-1 .and. j .le. ny/2-1 .and. k .le. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/(i**2+j**2+k**2)
             else if (i .le. nx/2-1 .and. j .le. ny/2-1 .and. k .gt. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/(i**2+j**2+(k-nz)**2)
             else if (i .le. nx/2-1 .and. j .gt. ny/2-1 .and. k .le. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/(i**2+(j-ny)**2+k**2)
             else if (i .gt. nx/2-1 .and. j .le. ny/2-1 .and. k .le. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/((i-nx)**2+j**2+k**2)
             else if (i .le. nx/2-1 .and. j .gt. ny/2-1 .and. k .gt. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/(i**2 + (j-ny)**2 + (k-nz)**2)
             else if (i .gt. nx/2-1 .and. j .gt. ny/2-1 .and. k .le. nz/2-1) then 
                output(i,j,k) = -input(i,j,k)/((i-nx)**2 + (j-ny)**2 + k**2)
             else if (i .gt. nx/2-1 .and. j .le. ny/2-1 .and. k .gt. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/((i-nx)**2 + j**2 + (k-nz)**2) 
             else if (i .gt. nx/2-1 .and. j .gt. ny/2-1 .and. k .gt. nz/2-1) then
                output(i,j,k) = -input(i,j,k)/((i-nx)**2 + (j-ny)**2 + (k-nz)**2)
             end if
          end do
       end do
    end do
    call x_derivative(output, nx, ny, nz, 0, 1, 0)
    call y_derivative(output, nx, ny, nz, 0, 1, 0)
    call z_derivative(output, nx, ny, nz, 0, 1, 0)
    call x_derivative(input, nx, ny, nz, 0, 1, 0)
    call y_derivative(input, nx, ny, nz, 0, 1, 0)
    call z_derivative(input, nx, ny, nz, 0, 1, 0)							
  end subroutine poisson

  subroutine double_deriv(input, nx, ny, nz, deriv1, deriv2)

    integer :: nx, ny, nz, deriv1, deriv2
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: input	

    if (deriv1 .eq. 1 .and. deriv2 .eq. 1) then
       call x_derivative(input, nx, ny, nz, 1, 1, 1)
       call x_derivative(input, nx, ny, nz, 1, 1, 1)
    elseif (deriv1 .eq. 1 .and. deriv2 .eq. 2 .or. deriv1 .eq. 2 .and. deriv2 .eq. 1) then
       call x_derivative(input, nx, ny, nz, 1, 1, 1)
       call y_derivative(input, nx, ny, nz, 1, 1, 1)
    elseif (deriv1 .eq. 1 .and. deriv2 .eq. 3 .or. deriv1 .eq. 3 .and. deriv2 .eq. 1) then
       call x_derivative(input, nx, ny, nz, 1, 1, 1)
       call z_derivative(input, nx, ny, nz, 1, 1, 1)
    elseif (deriv1 .eq. 2 .and. deriv2 .eq. 3 .or. deriv1 .eq. 3 .and. deriv2 .eq. 2) then	
       call y_derivative(input, nx, ny, nz, 1, 1, 1)
       call z_derivative(input, nx, ny, nz, 1, 1, 1)
    elseif (deriv1 .eq. 2 .and. deriv2 .eq. 2) then
       call y_derivative(input, nx, ny, nz, 1, 1, 1)
       call y_derivative(input, nx, ny, nz, 1, 1, 1)
    elseif (deriv1 .eq. 3 .and. deriv2 .eq. 3) then
       call z_derivative(input, nx, ny, nz, 1, 1, 1)
       call z_derivative(input, nx, ny, nz, 1, 1, 1)
    end if

  end subroutine double_deriv

  subroutine getrms(f, nx, ny, nz, fps)

    implicit none

    double complex f(0:nx-1,0:ny-1,0:nz-1)
    integer i,j,k, nx, ny, nz
    double precision fbar,fp,fps

    fbar=0.0;fp=0.0;fps=0.0

    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx-1
             fbar=fbar+real(f(i,j,k))
          end do
       end do
    end do

    fbar=fbar/(nx*ny*nz)

    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx-1
             fp=f(i,j,k)-fbar
             fps=fps+fp**2
          end do
       end do
    end do

    fps = fps/(nx*ny*nz)

    fps = dsqrt(fps)

  end subroutine getrms

  subroutine compressibility_calculation(uc, vc, wc, ui, vi, wi, nx, ny, nz, compressibility)

    integer :: nx, ny, nz
    double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: uc, vc ,wc, ui, vi, wi
    double precision :: compressibility, rmsval, compressible_sum, incompressible_sum

    compressible_sum = 0.0d0
    call getrms(uc, nx, ny, nz, rmsval)
    compressible_sum = compressible_sum + rmsval
    call getrms(vc, nx, ny, nz, rmsval)
    compressible_sum = compressible_sum + rmsval
    call getrms(wc, nx, ny, nz, rmsval)
    compressible_sum = compressible_sum + rmsval

    incompressible_sum = 0.0d0
    call getrms(ui, nx, ny, nz, rmsval)
    incompressible_sum = incompressible_sum + rmsval
    call getrms(vi, nx, ny, nz, rmsval)
    incompressible_sum = incompressible_sum + rmsval
    call getrms(wi, nx, ny, nz, rmsval)
    incompressible_sum = incompressible_sum + rmsval

    compressibility = compressible_sum/(incompressible_sum + compressible_sum)

  end subroutine compressibility_calculation

  subroutine wavenumbers_in_fourier_space(kx, ky, kz, mk, nx, ny, nz)

    double precision, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: mk
    double precision, dimension(0:nx-1) :: kx
    double precision, dimension(0:ny-1) :: ky
    double precision, dimension(0:nz-1) :: kz
    integer :: i, j, k
    integer :: nx, ny, nz

    do k = 0, nz-1
       kz(k) = k - nz/2
    end do

    do j = 0, ny-1
       ky(j) = j - ny/2
    end do

    do i = 0, nx-1
       kx(i) = i - nx/2
    end do

    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx-1
             mk(i,j,k) = dsqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
          end do
       end do
    end do

  end subroutine wavenumbers_in_fourier_space

  subroutine splitting_velocity_field(u, v, w, uc, vc, wc, ui, vi, wi, mk, kx, ky, kz, nx, ny, nz)

    integer :: nx, ny, nz, i, j, k
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: u, v, w
    double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: uc, vc, wc, ui, vi, wi
    double precision, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: mk
    double precision, dimension(0:nx-1) :: kx
    double precision, dimension(0:ny-1) :: ky
    double precision, dimension(0:nz-1) :: kz


    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx-1

             if ( mk(i,j,k) .ne. 0 ) then
                ui(i,j,k) = u(i,j,k) - ( kx(i)*u(i,j,k) + ky(j)*v(i,j,k) + &
                     & kz(k)*w(i,j,k) )*kx(i)/(mk(i,j,k))**2
                vi(i,j,k) = v(i,j,k) - ( kx(i)*u(i,j,k) + ky(j)*v(i,j,k) + &
                     & kz(k)*w(i,j,k) )*ky(j)/(mk(i,j,k))**2
                wi(i,j,k) = w(i,j,k) - ( kx(i)*u(i,j,k) + ky(j)*v(i,j,k) + &
                     & kz(k)*w(i,j,k) )*kz(k)/(mk(i,j,k))**2
             else if ( mk(i,j,k) .eq. 0 ) then
                ui(i,j,k) = 0.0d0
                vi(i,j,k) = 0.0d0
                wi(i,j,k) = 0.0d0
             end if

          end do
       end do
    end do

    uc = u - ui
    vc = v - vi
    wc = w - wi

  end subroutine splitting_velocity_field

end module ns

program compress

  use ns

  implicit none

  integer, parameter :: nx=96, ny=96, nz=96
  double precision, parameter :: pi=acos(-1.0d0), dx=2.0d0*pi/nx, dy=2.0d0*pi/ny, dz=2.0d0*pi/nz, gamma = 1.40d0, Mno = 0.30d0
  double precision, parameter :: tbase = 300.0d0, robase = 1.20d0, R = 287.150d0, pbase = robase*R*tbase, cbase = dsqrt(gamma*R*tbase)
  double precision, parameter :: t0 = 15000.0d0, ro0 = 1.20d0, p0 = ro0*R*t0, c0 = dsqrt(gamma*R*t0)
  double precision :: tempv1, tempv2, tempv3, mt
  double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: temp1, temp2, temp3, tempt, tempro, rhs_p, rhs_pt, p1, ro1, t1, pt, d, d_check
  ! double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: test, result
  double complex, dimension(0:nx-1,0:ny-1,0:nz-1) :: v1, v2, v3, w1, w2, w3, u1, u2, u3, u1_c, u2_c, u3_c, u1_i, u2_i, u3_i, du, dv, dw
  double precision :: dummy, v_rms, u_rms, q0, velocity_scaling, temperature_scaling, density_scaling, t_rms, p_rms, rho_rms, compressibility, divergence, idp
  integer :: i,j,k,statusint, temp_int, rms_dim_fluc
  double precision, dimension(0:nx-1) :: kx
  double precision, dimension(0:ny-1) :: ky
  double precision, dimension(0:nz-1) :: kz
  double precision, dimension(0:nx-1,0:ny-1,0:nz-1) :: mk
  double precision, dimension(0:83) :: e_u1, e_u2, e_u3, e_T, e_ro

  call wavenumbers_in_fourier_space(kx, ky, kz, mk, nx, ny, nz)

  statusint = 0

  rhs_p = 0.0d0
  temp1 = 0.0d0
  temp2 = 0.0d0
  v1 = 0.0d0
  v2 = 0.0d0
  v3 = 0.0d0

  open(unit=12,file="serial_initial.dat")
  do k=0,nz-1
     do j=0,ny-1
        do i=0,nx-1
           read(12,'(1x, 6E30.18)') dummy, dummy, dummy, tempv1, tempv2, tempv3
           v1(i,j,k)=tempv1
           v2(i,j,k)=tempv2
           v3(i,j,k)=tempv3
        end do
     end do
  end do
  close(12)

!  ! Dimensionalize the velocity fluctuations
!  v1 = v1*cbase
!  v2 = v2*cbase
!  v3 = v3*cbase

!  ! Re non-dimensionalize based on c0
!  v1 = v1/c0
!  v2 = v2/c0
!  v3 = v3/c0

  ! Calculating initial turbulent mach number
  q0 = 0.0d0
  call getrms(v1, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(v2, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(v3, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2

  mt = dsqrt(q0)

  write(*,*) "turbulent mach number before weakly compressible procedure is", mt

!  velocity_scaling = Mno/mt

!  v1 = v1*velocity_scaling
!  v2 = v2*velocity_scaling
!  v3 = v3*velocity_scaling

!  q0 = 0.0d0
!  call getrms(v1, nx, ny, nz, v_rms)
!  q0 = q0 + v_rms**2
!  call getrms(v2, nx, ny, nz, v_rms)
!  q0 = q0 + v_rms**2
!  call getrms(v3, nx, ny, nz, v_rms)
!  q0 = q0 + v_rms**2

!  write(*,*) "turbulent mach number before weakly compressible procedure after velocity_scaling is", dsqrt(q0)

  ! Dimensionalize velocity fluctuations
  v1 = v1*c0
  v2 = v2*c0
  v3 = v3*c0

  ! Calculating rms value for dimensional fluctuations
  q0 = 0.0d0
  call getrms(v1, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(v2, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(v3, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2

  rms_dim_fluc = dsqrt(q0)

  ! Non-dimensionalizing fluctuations based on rms dimensional fluctuations
  v1 = v1/rms_dim_fluc
  v2 = v2/rms_dim_fluc
  v3 = v3/rms_dim_fluc

  !call exit(statusint)

  temp1 = -1.0d0*v1*v1
  call double_deriv(temp1, nx, ny, nz, 1, 1)
  rhs_p = rhs_p + temp1

  temp1 = -1.0d0*v1*v2
  call double_deriv(temp1, nx, ny, nz, 1, 2)
  rhs_p = rhs_p + 2*temp1

  temp1 = -1.0d0*v1*v3
  call double_deriv(temp1, nx, ny, nz, 1, 3)
  rhs_p = rhs_p + 2*temp1

  temp1 = -1.0d0*v2*v2
  call double_deriv(temp1, nx, ny, nz, 2, 2)
  rhs_p = rhs_p + temp1

  temp1 = -1.0d0*v2*v3
  call double_deriv(temp1, nx, ny, nz, 2, 3)
  rhs_p = rhs_p + 2*temp1

  temp1 = -1.0d0*v3*v3
  call double_deriv(temp1, nx, ny, nz, 3, 3)
  rhs_p = rhs_p + temp1

  rhs_p = rhs_p

  call poisson(rhs_p, p1, nx, ny, nz)

  ! calculating density and temperature fluctuations from pressure fluctuations.
  ro1 = p1*(1.0d0/gamma)
  t1 = ((gamma-1.0d0)/(gamma))*p1

  rhs_pt = 0.0d0

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 1)
  temp2 = temp2 + temp1*v1
  temp1 = p1
  call x_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v1
  call double_deriv(temp2, nx, ny, nz, 1, 1)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 1)
  temp2 = temp2 + temp1*v2
  temp1 = p1
  call x_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v2
  call double_deriv(temp2, nx, ny, nz, 1, 2)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 1)
  temp2 = temp2 + temp1*v3
  temp1 = p1
  call x_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v3
  call double_deriv(temp2, nx, ny, nz, 1, 3)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 2)
  temp2 = temp2 + temp1*v1
  temp1 = p1
  call y_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v1
  call double_deriv(temp2, nx, ny, nz, 2, 1)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 2)
  temp2 = temp2 + temp1*v2
  temp1 = p1
  call y_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v2
  call double_deriv(temp2, nx, ny, nz, 2, 2)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 2)
  temp2 = temp2 + temp1*v3
  temp1 = p1
  call y_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v3
  call double_deriv(temp2, nx, ny, nz, 2, 3)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 3)
  temp2 = temp2 + temp1*v1
  temp1 = p1
  call z_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v1
  call double_deriv(temp2, nx, ny, nz, 3, 1)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 3)
  temp2 = temp2 + temp1*v2
  temp1 = p1
  call z_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v2
  call double_deriv(temp2, nx, ny, nz, 3, 2)
  rhs_pt = rhs_pt + temp2

  temp2 = 0.0d0
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 3)
  temp2 = temp2 + temp1*v3
  temp1 = p1
  call z_derivative(temp1, nx, ny, nz, 1, 1, 1)
  temp2 = temp2 + temp1*v3
  call double_deriv(temp2, nx, ny, nz, 3, 3)
  rhs_pt = rhs_pt + temp2

  rhs_pt = rhs_pt*2.0d0

  call poisson(rhs_pt, pt, nx, ny, nz)

  ! Calculating Dilatation d
  call vdotdel(nx, ny, nz, v1, v2, v3, p1, temp1, 4)
  d = (pt + temp1)*(-1.0d0/gamma)

  ! Calculating w
  temp1 = d
  call x_derivative(temp1, nx, ny, nz, 1, 0, 0)
  call y_derivative(temp1, nx, ny, nz, 1, 0, 0)
  call z_derivative(temp1, nx, ny, nz, 1, 0, 0)
  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1
           if (i .eq. 0 .and. j .eq. 0 .and. k .eq. 0) then
              w1(i,j,k) = 0.0
              w2(i,j,k) = 0.0
              w3(i,j,k) = 0.0
           else if (i .le. nx/2-1 .and. j .le. ny/2-1 .and. k .le. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-i/(i**2 + j**2 + k**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-j/(i**2 + j**2 + k**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-k/(i**2 + j**2 + k**2))*temp1(i,j,k)
           else if (i .le. nx/2-1 .and. j .le. ny/2-1 .and. k .gt. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-i/(i**2 + j**2 + (k-nz)**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-j/(i**2 + j**2 + (k-nz)**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-(k-nz)/(i**2 + j**2 + (k-nz)**2))*temp1(i,j,k)	
           else if (i .le. nx/2-1 .and. j .gt. ny/2-1 .and. k .le. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-i/(i**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-(j-ny)/(i**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-k/(i**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
           else if (i .gt. nx/2-1 .and. j .le. ny/2-1 .and. k .le. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-(i-nx)/((i-nx)**2 + j**2 + k**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-j/((i-nx)**2 + j**2 + k**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-k/((i-nx)**2 + j**2 + k**2))*temp1(i,j,k)
           else if (i .le. nx/2-1 .and. j .gt. ny/2-1 .and. k .gt. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-i/(i**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-(j-ny)/(i**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-(k-nz)/(i**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
           else if (i .gt. nx/2-1 .and. j .gt. ny/2-1 .and. k .le. nz/2-1) then 
              w1(i,j,k) = dcmplx(0.0d0,-(i-nx)/((i-nx)**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-(j-ny)/((i-nx)**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-k/((i-nx)**2 + (j-ny)**2 + k**2))*temp1(i,j,k)
           else if (i .gt. nx/2-1 .and. j .le. ny/2-1 .and. k .gt. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-(i-nx)/((i-nx)**2 + j**2 + (k-nz)**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-j/((i-nx)**2 + j**2 + (k-nz)**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-(k-nz)/((i-nx)**2 + j**2 + (k-nz)**2))*temp1(i,j,k)
           else if (i .gt. nx/2-1 .and. j .gt. ny/2-1 .and. k .gt. nz/2-1) then
              w1(i,j,k) = dcmplx(0.0d0,-(i-nx)/((i-nx)**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
              w2(i,j,k) = dcmplx(0.0d0,-(j-ny)/((i-nx)**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
              w3(i,j,k) = dcmplx(0.0d0,-(k-nz)/((i-nx)**2 + (j-ny)**2 + (k-nz)**2))*temp1(i,j,k)
           end if
        end do
     end do
  end do

  call x_derivative(w1, nx, ny, nz, 0, 1, 0)
  call y_derivative(w1, nx, ny, nz, 0, 1, 0)
  call z_derivative(w1, nx, ny, nz, 0, 1, 0)	

  call x_derivative(w2, nx, ny, nz, 0, 1, 0)
  call y_derivative(w2, nx, ny, nz, 0, 1, 0)
  call z_derivative(w2, nx, ny, nz, 0, 1, 0)	

  call x_derivative(w3, nx, ny, nz, 0, 1, 0)
  call y_derivative(w3, nx, ny, nz, 0, 1, 0)
  call z_derivative(w3, nx, ny, nz, 0, 1, 0)	

  ! write u and w to file
  u1 = v1 + gamma*Mno**2*w1
  u2 = v2 + gamma*Mno**2*w2
  u3 = v3 + gamma*Mno**2*w3

!  call splitting_velocity_field(u1, u2, u3, u1_c, u2_c, u3_c, u1_i, u2_i, u3_i, mk, kx, ky, kz, nx, ny, nz)
  u1_c = gamma*Mno**2*w1
  u2_c = gamma*Mno**2*w2
  u3_c = gamma*Mno**2*w3

  u1_i = v1
  u2_i = v2
  u3_i = v3

  ! Calculating dilatation of incompressibility
  du = u1_i
  call x_derivative(du, nx, ny, nz, 1, 1, 1)
  dv = u2_i
  call y_derivative(dv, nx, ny, nz, 1, 1, 1)
  dw = u3_i
  call z_derivative(dw, nx, ny, nz, 1, 1, 1)

  divergence = sum((real(du)+real(dv)+real(dw))**2)/(1.0d0*nx*ny*nz)
  write(*,*) "dilatation of incompressible part of flow-field is", divergence

  call compressibility_calculation(u1_c, u2_c, u3_c, u1_i, u2_i, u3_i, nx, ny, nz, compressibility)
  write(*,*) "compressibility in flow field is", compressibility

  !call exit(statusint)

  u1 = u1*rms_dim_fluc
  u2 = u2*rms_dim_fluc
  u3 = u3*rms_dim_fluc

  u1 = u1/c0
  u2 = u2/c0
  u3 = u3/c0

  ! re- scaling p and t.
  t1 = t1*gamma*Mno**2
  call getrms(t1, nx, ny, nz, t_rms)
  write(*,*) "t_rms", t_rms

  p1 = p1*gamma*Mno**2
  p1 = p1*p0/(ro0*c0**2)
  call getrms(p1, nx, ny, nz, p_rms)
  write(*,*) "p_rms", p_rms

  ro1 = ro1*gamma*Mno**2
  call getrms(ro1, nx, ny, nz, rho_rms)
  write(*,*) "rho_rms", rho_rms

  q0 = 0.0d0
  call getrms(u1, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(u2, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2
  call getrms(u3, nx, ny, nz, v_rms)
  q0 = q0 + v_rms**2

  write(*,*) "turbulent mach number after weakly compressible flow field is", dsqrt(q0)

  e_u1 = 0.0d0
  e_u2 = 0.0d0
  e_u3 = 0.0d0
  e_ro = 0.0d0
  e_T = 0.0d0

  temp1 = u1
  call x_derivative(temp1, nx, ny, nz, 1, 0, 0)
  call y_derivative(temp1, nx, ny, nz, 1, 0, 0)
  call z_derivative(temp1, nx, ny, nz, 1, 0, 0)

  temp2 = u2
  call x_derivative(temp2, nx, ny, nz, 1, 0, 0)
  call y_derivative(temp2, nx, ny, nz, 1, 0, 0)
  call z_derivative(temp2, nx, ny, nz, 1, 0, 0)

  temp3 = u3
  call x_derivative(temp3, nx, ny, nz, 1, 0, 0)
  call y_derivative(temp3, nx, ny, nz, 1, 0, 0)
  call z_derivative(temp3, nx, ny, nz, 1, 0, 0)

  tempro = ro1
  call x_derivative(tempro, nx, ny, nz, 1, 0, 0)
  call y_derivative(tempro, nx, ny, nz, 1, 0, 0)
  call z_derivative(tempro, nx, ny, nz, 1, 0, 0)

  tempt = t1
  call x_derivative(tempt, nx, ny, nz, 1, 0, 0)
  call y_derivative(tempt, nx, ny, nz, 1, 0, 0)
  call z_derivative(tempt, nx, ny, nz, 1, 0, 0)

  call wrap_around_3d(temp1, nx, ny, nz)
  call wrap_around_3d(temp2, nx, ny, nz)
  call wrap_around_3d(temp3, nx, ny, nz)
  call wrap_around_3d(tempro, nx, ny, nz)
  call wrap_around_3d(tempt, nx, ny, nz)

  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1

           if ( mk(i,j,k) .ne. 0 ) then

              if ( mk(i,j,k) .gt. int(mk(i,j,k))+0.50d0 .and. mk(i,j,k) .lt. nint(mk(i,j,k))+0.50d0 ) then
                 temp_int = nint(mk(i,j,k))           
              else
                 temp_int = int(mk(i,j,k))
              end if

              e_u1(temp_int) = e_u1(temp_int) + temp1(i,j,k)*dconjg(temp1(i,j,k))
              e_u2(temp_int) = e_u2(temp_int) + temp2(i,j,k)*dconjg(temp2(i,j,k))
              e_u3(temp_int) = e_u3(temp_int) + temp3(i,j,k)*dconjg(temp3(i,j,k))
              e_ro(temp_int) = e_ro(temp_int) + tempro(i,j,k)*dconjg(tempro(i,j,k))
              e_T(temp_int) = e_T(temp_int) + tempt(i,j,k)*dconjg(tempt(i,j,k))

           else

              e_u1(temp_int) = 0.0d0
              e_u2(temp_int) = 0.0d0
              e_u3(temp_int) = 0.0d0
              e_ro(temp_int) = 0.0d0
              e_T(temp_int) = 0.0d0

           end if

        end do
     end do
  end do

  open(unit=54, file="spectra_velocity_ristorcelli.dat")
  write(54,*) 'variables = "k", "e_u1", "e_u2", "e_u3","e_ro","e_T"'
  do i = 0, 83
     idp = 1.0d0*i
     write(54,'(6E30.18)') idp, e_u1(i), e_u2(i), e_u3(i), e_ro(i), e_T(i)
  end do
  close(54)

  open(unit=23,file="serial_initial2.dat")
  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1
           write(23,'(6E30.18)') real(u1(i,j,k)), real(u2(i,j,k)), real(u3(i,j,k)), real(p1(i,j,k)), real(ro1(i,j,k)), real(t1(i,j,k))
        end do
     end do
  end do
  close(23)

end program compress

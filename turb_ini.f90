module ns

	implicit none

	include "fftw3.f"

	contains

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

	end subroutine


	subroutine x_derivative(input_3d, nx, ny, nz)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:nx-1) :: input_x
		integer*8 ::plan_x_fwd, plan_x_bwd

		call dfftw_plan_dft_1d(plan_x_fwd, nx, input_x, input_x, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_x_bwd, nx, input_x, input_x, fftw_backward, fftw_measure)
		

		do k = 0, nz-1
			do j = 0, ny-1										! transforming to fourier space.
				input_x = input_3d(:, j, k)							! copy to 1d array from 3d array.
				call dfftw_execute_dft(plan_x_fwd, input_x, input_x)				! call fftw 1d fft subroutine.
				call fourier_derivative(nx, input_x)						! call derivative subroutine.
				input_3d(:, j, k) = input_x							! copy back to 3d input array from
			end do											! 1d array.
		end do

		input_3d = input_3d/nx										! renormalization after x-fft.

		do k = 0, nz-1
			do j = 0, ny-1										! transforming back to physical space.
				input_x = input_3d(:, j, k)
				call dfftw_execute_dft(plan_x_bwd, input_x, input_x)
				input_3d(:, j, k) = input_x
			end do
		end do

		call dfftw_destroy_plan(plan_x_fwd)
		call dfftw_destroy_plan(plan_x_bwd)

	end subroutine


	subroutine y_derivative(input_3d, nx, ny, nz)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:ny-1) :: input_y
		integer*8 ::plan_y_fwd, plan_y_bwd

		call dfftw_plan_dft_1d(plan_y_fwd, ny, input_y, input_y, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_y_bwd, ny, input_y, input_y, fftw_backward, fftw_measure)
		

		do k = 0, nz-1
			do i = 0, nx-1
				input_y = input_3d(i, :, k)
				call dfftw_execute_dft(plan_y_fwd, input_y, input_y)
				call fourier_derivative(ny, input_y)
				input_3d(i, :, k) = input_y
			end do
		end do

		input_3d = input_3d/ny

		do k = 0, nz-1
			do i = 0, nx-1
				input_y = input_3d(i, :, k)
				call dfftw_execute_dft(plan_y_bwd, input_y, input_y)
				input_3d(i, :, k) = input_y
			end do
		end do

		call dfftw_destroy_plan(plan_y_fwd)
		call dfftw_destroy_plan(plan_y_bwd)

	end subroutine


	subroutine z_derivative(input_3d, nx, ny, nz)
		
		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:nz-1) :: input_z
		integer*8 ::plan_z_fwd, plan_z_bwd

		call dfftw_plan_dft_1d(plan_z_fwd, nz, input_z, input_z, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_z_bwd, nz, input_z, input_z, fftw_backward, fftw_measure)
		

		do j = 0, ny-1
			do i = 0, nx-1
				input_z = input_3d(i, j, :)
				call dfftw_execute_dft(plan_z_fwd, input_z, input_z)
				call fourier_derivative(nz, input_z)
				input_3d(i, j, :) = input_z
			end do
		end do

		input_3d = input_3d/nz

		do j = 0, ny-1
			do i = 0, nx-1
				input_z = input_3d(i, j, :)
				call dfftw_execute_dft(plan_z_bwd, input_z, input_z)
				input_3d(i, j, :) = input_z
			end do
		end do

		call dfftw_destroy_plan(plan_z_fwd)
		call dfftw_destroy_plan(plan_z_bwd)

	end subroutine

end module

program parallel

	use ns

	implicit none

	integer :: i, j, k, n, nnpx, p
	integer, parameter :: nx=96, ny=96, nz=96
	integer, parameter :: nproc = 4, nproot = 2, npx = nx/nproot
	double precision, parameter :: Re = 515.0d0
	character (len=50) :: filename
	character (len=11), parameter :: const = "Turb_inflow"
	character (len=4), parameter :: dat = ".dat"
	character (len=1024) :: dummy
	double precision :: xt,yt,zt,ut,vt,wt,u_avg, v_avg, w_avg,q, lambx, lamby, lambz, lamb, Rel, sx, sy, sz, skew
	double precision, dimension(0:nx-1,0:ny-1,0:nz-1) :: u, v, w
	double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: temp, div

	div = (0, 0)

	do n = 0, nproc-1
		if (n .lt. 10) then
			write(filename, '(a,i0,i0,a)') const,0,n,dat
			open(unit=n+10,file=trim(filename))
		else
			write(filename, '(a,i0,a)') const,n,dat
			open(unit=n+10,file=trim(filename))
		end if
	end do

	open(unit=9,file="serial_initial.dat")
	do n = 0, nproc-1
		read(n+10,*) dummy
		read(n+10,*) dummy
	end do
	do p = 0, nproot-1
		do nnpx = 1, npx
			do n = p, nproc-1, nproot
				i = 0
				do while (i .le. (ny*npx)-1)
					read(n+10,'(1x,6E30.18)') xt, yt, zt, ut, vt, wt
					write(9,'(1x,6E30.18)') xt, yt, zt, ut, vt, wt
					i = i+1
				end do
			end do
		end do
	end do
	do n = 0, nproc-1
		close(n+10)
	end do
	close(9)

	open(unit=9,file="serial_initial.dat")
	do k = 0, nz-1
		do j = 0, ny-1	
			do i = 0, nx-1
				read(9, '(1x,6E30.18)') xt, yt, zt, u(i,j,k), v(i,j,k), w(i,j,k)
			end do
		end do
	end do
	close(9)	

	u_avg = sum(u)/(nx*ny*nz)
	v_avg = sum(v)/(nx*ny*nz)
	w_avg = sum(w)/(nx*ny*nz)	

	print *, sum(u)/(nx*ny*nz), sum(v)/(nx*ny*nz), sum(w)/(nx*ny*nz)

	print *, sum(u**2)/(nx*ny*nz), sum(v**2)/(nx*ny*nz), sum(w**2)/(nx*ny*nz)
	
	temp = u
	call x_derivative(temp, nx, ny, nz)
	div = div + temp
	lambx = dsqrt(sum(u**2)/sum(real(temp)**2))
	sx = sum(real(temp)**3)/(nx*ny*nz)/(sum(real(temp)**2)/(nx*ny*nz))**1.50d0

	temp = v
	call y_derivative(temp, nx, ny, nz)
	div = div + temp
	lamby = dsqrt(sum(v**2)/sum(real(temp)**2))
	sy = sum(real(temp)**3)/(nx*ny*nz)/(sum(real(temp)**2)/(nx*ny*nz))**1.50d0
	
	temp = w
	call z_derivative(temp, nx, ny, nz)
	div = div + temp
	lambz = dsqrt(sum(w**2)/sum(real(temp)**2))
	sz = sum(real(temp)**3)/(nx*ny*nz)/(sum(real(temp)**2)/(nx*ny*nz))**1.50d0

	lamb = (lambx+lamby+lambz)/3.0d0
	skew = (sx+sy+sz)/3.0d0

	q = (sum(u**2)+sum(v**2)+sum(w**2))/(nx*ny*nz)

	Rel = Re*dsqrt(q)*lamb/dsqrt(3.0d0)

	print *, "q1 value is", sqrt(q)
	print *, "averaged taylor microscale is", lamb
	print *, "Taylor Reynolds Number is", Rel
	print *, "Max Velocity Divergence is", maxval(abs(real(div)))
	print *, "Velocity Derivative Skewness is", skew

	

end program



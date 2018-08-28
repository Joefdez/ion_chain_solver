
! module containing variable initializations

module initialization
  use constants
  implicit none

contains
  ! based on http://jblevins.org/log/control-file
  subroutine initialize_system(n_particles, mass1, charge1, tt, pt, dt, traj, nsave, &
             long_freq, alphay, alphaz, ic_radius, savefreq)!, initT)

    implicit none

    character(len=100) :: buffer!, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0

    integer, intent(inout)            :: n_particles   ! dimensionality of the problem, number of partilcles
    real(kind=8), intent(inout)       :: mass1, charge1, ic_radius
    real(kind=8), intent(inout)       :: alphay, alphaz, long_freq!, initT
    real(kind=8), intent(inout)       :: dt, tt ,pt
    integer, intent(inout)            :: traj
    integer, intent(inout)            :: nsave, savefreq


    open(unit=fh, file='values.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first isntance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
          case('n_particles')
              print*, "n_particles"
              read(buffer,*,iostat=ios) n_particles
              print*, "nparticles", n_particles
          case('mass1')
              read(buffer,*,iostat=ios) mass1
              print*, "mass1", mass1
          case('charge1')
              read(buffer,*,iostat=ios) charge1
              print*, "charge1", charge1
          case('dt')
              read(buffer,*,iostat=ios) dt
              print*, "dt", dt
          case('tt')
              read(buffer,*,iostat=ios) tt
              print*, "tt", tt
          case('pt')
              read(buffer,*,iostat=ios) pt
              print*, "pt", pt
          case('traj')
              read(buffer,*,iostat=ios) traj
              print*, "traj", traj
          case('alphay')
              read(buffer,*,iostat=ios) alphay
              print*, "alphay", alphay
          case('alphaz')
                read(buffer,*,iostat=ios) alphaz
                print*, "alphaz", alphaz
          case('icradius')
              read(buffer,*,iostat=ios) ic_radius
              print*, "icradius", ic_radius
          case('long_freq')
              read(buffer,*,iostat=ios) long_freq
              print*, "long_freq", long_freq
          case('nsave')
              read(buffer,*,iostat=ios) nsave
              print*, "nsave", nsave
          case('savefreq')
              read(buffer,*,iostat=ios) savefreq
              print*, "save_freq", savefreq
          case default
              print*, "Skipping invalid value."
        end select
      end if
    end do

    print*, nsave

    close(unit=fh)

  end subroutine initialize_system

  subroutine initialize_laser_chain(delx1, dely1, delz1, delx2, dely2, delz2, Gam, omega, &
                                    Ix1, Iy1, Iz1, Ix2, Iy2, Iz2)
    ! initialize laser and ion interaction parameters

    implicit none

    character(len=100) :: buffer!, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0
    real(kind=8), intent(inout) :: delx1, dely1, delz1, delx2, dely2, delz2
    real(kind=8), intent(inout)  :: Gam, omega
    real(kind=8), intent(inout) :: Ix1, Iy1, Iz1, Ix2, Iy2, Iz2

    open(unit=fh, file='chain_laser.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first isntance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
        case('detuning_left')
              read(buffer,*,iostat=ios) delx1, dely1, delz1
          case('detuning_right')
              read(buffer,*,iostat=ios) delx2, dely2, delz2
          case('linewidth')
              read(buffer,*,iostat=ios) Gam
          case('target_frequency')
              read(buffer,*,iostat=ios) omega
          case('intensity_left')
              read(buffer,*,iostat=ios) Ix1, Iy1, Iz1
          case('intensity_right')
              read(buffer,*,iostat=ios) Ix2, Iy2, Iz2
          case default
              print*, "Skipping invalid value."
        end select
      end if
    end do

    print*,"read"
  end subroutine initialize_laser_chain

  subroutine dimensionless_doppler_values(eta, D, mass, long_freq, char_length, aeta, aD)
    implicit none
    real(kind=8), intent(in)              :: eta, D, mass, long_freq, char_length
    real(kind=8), intent(inout)           :: aeta, aD
    ! calculate scaled friction and difussion coefficients
    aeta = eta / (mass * long_freq)
    aD   = D / (char_length*char_length*mass*mass*long_freq*long_freq*long_freq)

  end subroutine dimensionless_doppler_values


  subroutine initcond(qq, pp)
    implicit none
    real(kind=8), dimension(:), intent(inout) :: qq, pp
    integer, parameter :: fh=15

    open(unit=fh,file='initcond.dat', action='read')
    read(fh,*) qq
    read(fh,*) pp


    close(unit=1)
  end subroutine initcond


  subroutine initLaser(mass, char_length, long_freq, aetax1, aetax2, aetay1, aetay2, aetaz1, aetaz2, &
                       aDx1, aDx2, aDy1, aDy2, aDz1, aDz2)
    implicit none
    real(kind=8), intent(in)    :: mass, char_length, long_freq
    real(kind=8), intent(inout) :: aetax1, aetax2, aetay1, aetay2, aetaz1, aetaz2
    real(kind=8), intent(inout) :: aDx1, aDx2, aDy1, aDy2, aDz1, aDz2
    real(kind=8) :: delx1, delx2, dely1, dely2, delz1, delz2, Gam, omega0, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2
    real(kind=8) :: etax1, etax2, etay1, etay2, etaz1, etaz2
    real(kind=8) :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    real(kind=8) :: kx1, kx2, ky1, ky2, kz1, kz2
    real(kind=8) :: omegax1, omegax2, omegay1, omegay2, omegaz1, omegaz2

    call initialize_laser_chain(delx1, dely1, delz1, delx2, dely2, delz2, Gam, omega0, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2)
    Gam  = 2.0d0 * pi * Gam         ! Natural linewidth
    delx1  = delx1 * Gam
    dely1  = dely1 * Gam
    delz1  = delz1 * Gam
    delx2  = delx2 * Gam
    dely2  = dely2 * Gam
    delz2  = delz2 * Gam
    omega0 = 2.0d0 * pi * omega0    ! Central frequency
    omegax1 = delx1 + omega0
    omegay1 = dely1 + omega0
    omegaz1 = delz1 + omega0
    omegax2 = delx2 + omega0
    omegay2 = dely2 + omega0
    omegaz2 = delz2 + omega0

    kx1 =  omegax1 / cc   ! Not really wavelengths, dimensionally they are wavenumbers
    ky1 =  omegay1 / cc
    kz1 =  omegaz1 / cc
    kx2 =  omegax2 / cc
    ky2 =  omegay2 / cc
    kz2 =  omegaz2 / cc
   ! Calculate diffusion and friction coefficients
    etax1 = -4.0d0*hbar*kx1*kx1*Ix1*(2.0d0*delx1/Gam)/( (1 + 4.0d0*delx1*delx1/(Gam*Gam)) * (1 + 4.0d0*delx1*delx1/(Gam*Gam)) )
    etay1 = -4.0d0*hbar*ky1*ky1*Iy1*(2.0d0*dely1/Gam)/( (1 + 4.0d0*dely1*dely1/(Gam*Gam)) * (1 + 4.0d0*dely1*dely1/(Gam*Gam)) )
    etaz1 = -4.0d0*hbar*kz1*kz1*Iz1*(2.0d0*delz1/Gam)/( (1 + 4.0d0*delz1*delz1/(Gam*Gam)) * (1 + 4.0d0*delz1*delz1/(Gam*Gam)) )
    etax1 = -4.0d0*hbar*kx2*kx2*Ix2*(2.0d0*delx2/Gam)/( (1 + 4.0d0*delx2*delx2/(Gam*Gam)) * (1 + 4.0d0*delx2*delx2/(Gam*Gam)) )
    etay2 = -4.0d0*hbar*ky2*ky2*Iy2*(2.0d0*dely2/Gam)/( (1 + 4.0d0*dely2*dely2/(Gam*Gam)) * (1 + 4.0d0*dely2*dely2/(Gam*Gam)) )
    etaz2 = -4.0d0*hbar*kz2*kz2*Iz2*(2.0d0*delz2/Gam)/( (1 + 4.0d0*delz2*delz2/(Gam*Gam)) * (1 + 4.0d0*delz2*delz2/(Gam*Gam)) )
    Dx1   = hbar*hbar*kx1*kx1*Ix1*(Gam)/(1.0d0 + 4.0d0*delx1*delx1/(Gam*Gam))
    Dy1   = hbar*hbar*ky1*ky1*Iy1*(Gam)/(1.0d0 + 4.0d0*dely1*dely1/(Gam*Gam))
    Dz1   = hbar*hbar*kz1*kz1*Iz1*(Gam)/(1.0d0 + 4.0d0*delz1*delz1/(Gam*Gam))
    Dx2   = hbar*hbar*kx2*kx2*Ix2*(Gam)/(1.0d0 + 4.0d0*delx2*delx2/(Gam*Gam))
    Dy2   = hbar*hbar*ky2*ky2*Iy2*(Gam)/(1.0d0 + 4.0d0*dely2*dely2/(Gam*Gam))
    Dz2   = hbar*hbar*kz2*kz2*Iz2*(Gam)/(1.0d0 + 4.0d0*delz2*delz2/(Gam*Gam))
    !initT = initT * (kb/(char_length*char_length*mass*long_freq*long_freq))
    !initSpeed = sqrt(2*initSpeed)

    ! Calculate dimensionless Doppler cooling parameters
    call dimensionless_doppler_values(etax1, Dx1, mass, long_freq, char_length, aetax1, aDx1)
    call dimensionless_doppler_values(etay1, Dy1, mass, long_freq, char_length, aetay1, aDy1)
    call dimensionless_doppler_values(etaz1, Dz1, mass, long_freq, char_length, aetaz1, aDz1)
    call dimensionless_doppler_values(etax2, Dx2, mass, long_freq, char_length, aetax2, aDx2)
    call dimensionless_doppler_values(etay2, Dy2, mass, long_freq, char_length, aetay2, aDy2)
    call dimensionless_doppler_values(etaz2, Dz2, mass, long_freq, char_length, aetaz2, aDz2)

  end subroutine initLaser

end module initialization

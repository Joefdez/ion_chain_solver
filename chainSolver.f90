program chainSolver

  use mpi
  use constants
  use randomModule
  use initcond_generator
  use initialization
  use support_functions_threed
  implicit none

  integer :: nparticles, traj, nsteps, nsave, ii, jj, ll, mm, saved, kk, nbath, save_freq, REM, ps, local_traj, ntrajC
  real(kind=8) :: mass, charge, tt, tp, dt, dst, long_freq, alphax, alphay, alphaz, ic_radius, initT, pt
  real(kind=8) :: aetax1, aetax2, aetay1, aetay2, aetaz1, aetaz2
  real(kind=8) :: aDx1, aDx2, aDy1, aDy2, aDz1, aDz2
  real(kind=8) :: JJix, JJiy, JJiz, JJix_s, JJiy_s, JJiz_s, JJix_av, JJiy_av, JJiz_av
  real(kind=8) :: char_length!, energy
  real(kind=8), dimension(:), allocatable  :: energy
  real(kind=8), dimension(:,:), allocatable :: hc!, energy
  real(kind=8), dimension(:), allocatable   :: Axx, Ayy, Azz, Axxi, Ayyi, Azzi
  real(kind=8), dimension(:), allocatable   :: Apx, Apy, Apz, Apxi, Apyi, Apzi
  real(kind=8), dimension(:), allocatable   :: stermsBx, stermsBy, stermsBz
  real(kind=8), dimension(:,:), allocatable :: fx1, fx2, fy1, fy2, fz1, fz2
  real(kind=8), dimension(:), allocatable   :: fx, fy, fz
  real(kind=8), dimension(:), allocatable   :: dOmx, dOmy, dOmz
  real(kind=8), dimension(:), allocatable   :: xx0, yy0, zz0
  real(kind=8), dimension(:), allocatable   :: xx, yy, zz, xxi, yyi, zzi, xxN, yyN, zzN
  real(kind=8), dimension(:), allocatable   :: px, py, pz, pxi, pyi, pzi, pxN, pyN, pzN
  real(kind=8), dimension(:,:), allocatable :: px2S, py2S, pz2S, px2S_av, py2S_av, pz2S_av
  real(kind=8), dimension(:,:), allocatable :: xxS, yyS, zzS, pxS, pyS, pzS
  real(kind=8), dimension(:,:), allocatable :: xxS_av, yyS_av, zzS_av, pxS_av, pyS_av, pzS_av
  real(kind=8), dimension(:,:), allocatable :: invD1,invD2
  logical                                   :: continue, opened
  ! mpi variables
  integer :: rank, procs, status(MPI_STATUS_SIZE), alloc_err, source, ierr


  call MPI_INIT(ierr)                                                                               ! Neccesary mpi initialization calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

  if(rank .eq. 0 ) then
    print*, "Reading solution and system parameters"
    call initialize_system(nparticles, mass, charge, tt, pt, dt, traj, nsave, &
                           long_freq, alphay, alphaz, ic_radius, save_freq)!, initT)
  end if

  call mpi_bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(traj, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nsteps, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(save_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nsave, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(mass, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(charge, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(tt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(pt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alphax, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alphay, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alphaz, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ic_radius, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(long_freq, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(initT, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(mpi_comm_world, ierr)

  mass   = mass*uu
  charge = charge*ee
  long_freq = 2.0d0*pi*long_freq
  dst = sqrt(dt)
  nsteps=int(tt/dt)
  char_length = ((charge*charge/(4.0d0*pi*ep0))/(mass*long_freq*long_freq))**(1.0/3.0)
  !nssteps = int(nsteps/save_freq)  ! Not saving every single timestep saves memory. Must ask about this
  !fin = 0.8d0*nsteps

  if(rank .eq. 0) then
    ! Initialize laser parameters
    ! The main code has access to the adimensional Doppler parameters
    print*, "Reading laser parameters"
    call initLaser(mass, char_length, long_freq, aetax1, aetax2, aetay1, aetay2, aetaz1, aetaz2, &
                     aDx1, aDx2, aDy1, aDy2, aDz1, aDz2)
  end if

  ! Share the Doppler parameters with the rest of the processes
  call mpi_bcast(aetax1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aetay1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aetaz1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aetax2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aetay2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aetaz2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDx1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDy1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDz1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDx2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDy2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(aDz2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1
  nbath = 3
  print*, nsave
  ! Allocate the different arrays
  allocate(xx(1:nparticles))
  xx = 0.0d0
  allocate(xxN(1:nparticles))
  xxN = 0.0d0
  allocate(yy(1:nparticles))
  yy = 0.0d0
  allocate(yyN(1:nparticles))
  yyN = 0.0d0
  allocate(zz(1:nparticles))
  zz = 0.0d0
  allocate(zzN(1:nparticles))
  zzN = 0.0d0
  allocate(xx0(1:nparticles))
  allocate(yy0(1:nparticles))
  allocate(zz0(1:nparticles))
  allocate(px(1:nparticles))
  allocate(pxN(1:nparticles))
  allocate(py(1:nparticles))
  allocate(pyN(1:nparticles))
  allocate(pz(1:nparticles))
  allocate(pzN(1:nparticles))
  allocate(xxS(1:nparticles, 1:nsave))
  allocate(yyS(1:nparticles, 1:nsave))
  allocate(zzS(1:nparticles, 1:nsave))
  allocate(pxS(1:nparticles, 1:nsave))
  allocate(pyS(1:nparticles, 1:nsave))
  allocate(pzS(1:nparticles, 1:nsave))
  allocate(xxS_av(1:nparticles, 1:nsave))
  allocate(yyS_av(1:nparticles, 1:nsave))
  allocate(zzS_av(1:nparticles, 1:nsave))
  allocate(px2s(1:nparticles, 1:nsave))
  px2s =0.0d0
  allocate(py2s(1:nparticles, 1:nsave))
  py2s = 0.0d0
  allocate(pz2s(1:nparticles, 1:nsave))
  pz2s = 0.0d0
  allocate(px2s_av(1:nparticles, 1:nsave))
  px2s_av =0.0d0
  allocate(py2s_av(1:nparticles, 1:nsave))
  py2s_av = 0.0d0
  allocate(pz2s_av(1:nparticles, 1:nsave))
  pz2s_av = 0.0d0

  ! Arrays for storing Coulmb forces
  allocate(fx1(1:nparticles,1:nparticles))
  fx1 = 0.0d0
  allocate(fx2(1:nparticles,1:nparticles))
  fx2 = 0.0d0
  allocate(fy1(1:nparticles,1:nparticles))
  fy1 = 0.0d0
  allocate(fy2(1:nparticles,1:nparticles))
  fy2 = 0.0d0
  allocate(fz1(1:nparticles,1:nparticles))
  fz1 = 0.0d0
  allocate(fz2(1:nparticles,1:nparticles))
  fz2 = 0.0d0
  allocate(invD1(1:nparticles,1:nparticles))
  invD1 = 0.0d0
  allocate(invD2(1:nparticles,1:nparticles))
  invD2 = 0.0d0
  allocate(fx(1:nparticles))
  fx = 0.0d0
  allocate(fy(1:nparticles))
  fy = 0.0d0
  allocate(fz(1:nparticles))
  fz = 0.0d0
  allocate(Axx(1:nparticles))
  ! Arrays for stochastic elements
  Axx = 0.0d0
  allocate(Ayy(1:nparticles))
  Ayy = 0.0d0
  allocate(Azz(1:nparticles))
  Azz = 0.0d0
  allocate(Apx(1:nparticles))
  Apx = 0.0d0
  allocate(Apy(1:nparticles))
  Apy = 0.0d0
  allocate(Apz(1:nparticles))
  Apz = 0.0d0
  allocate(Axxi(1:nparticles))
  Axxi = 0.0d0
  allocate(Ayyi(1:nparticles))
  Ayyi = 0.0d0
  allocate(Azzi(1:nparticles))
  Azzi = 0.0d0
  allocate(Apxi(1:nparticles))
  Apxi = 0.0d0
  allocate(Apyi(1:nparticles))
  Apyi = 0.0d0
  allocate(Apzi(1:nparticles))
  Apzi = 0.0d0
  allocate(dOmx(1:nparticles))
  dOmx = 0.0d0
  allocate(dOmy(1:nparticles))
  dOmy = 0.0d0
  allocate(dOmz(1:nparticles))
  dOmz = 0.0d0
  allocate(stermsBx(1:nparticles))
  stermsBx(1:nbath) = sqrt(2.0d0*aDx1)
  stermsBx((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aDx2)
  allocate(stermsBy(1:nparticles))
  stermsBy(1:nbath) = sqrt(2.0d0*aDy1)
  stermsBy((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aDy2)
  allocate(stermsBz(1:nparticles))
  stermsBz(1:nbath) = sqrt(2.0d0*aDz1)
  stermsBz((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aDz2)


  allocate(energy(1:nparticles))
  energy = 0.0d0
  allocate(hc(1:nparticles,1:nparticles))
  energy = 0.0d0
  print*, "allocations completed"

  ! Initial conditions
  if (rank .eq. 0) then
    open(unit=15, file="eq_pos.dat", action='read')
    do ii = 1, nparticles, 1
      read(15,*) xx0(ii), yy0(ii), zz0(ii)
    end do
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_bcast(xx0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(yy0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(zz0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  ! Solution of the equations of motion
  do kk=1, local_traj, 1
    print*, "Proc.", rank, "on trajectory", kk
    !xxs = 0.0d0
    !yys = 0.0d0
    !zzs = 0.0d0
    px = 0.0d0
    py = 0.0d0
    pz = 0.0d0
    ll = 0
    mm = 1
    JJix = 0.0d0
    JJiy = 0.0d0
    JJiz = 0.0d0
    saved = 1
    tp    = 0.0d0
    call icpgen(nparticles, 0.005d0, xx0, yy0, zz0, xx, yy, zz)
    !call icmomgen(nparticles, initSpeed, pxold, pyold, pzold)
    call ranseed()
!   JJix_av  = 0.0d0
!   JJiy_av = 0.0d0
    continue = .True.
    ps = 0
    do while(continue)
      !print*, tp, tt
      call coulombM(nparticles, xx, yy, zz, fx1, fy1, fz1, invD1)
      fx = 0.0d0
      fy = 0.0d0
      fz = 0.0d0
      do jj=1, nparticles, 1
        fx(jj) = sum(fx1(jj,:), 1)
        fy(jj) = sum(fy1(jj,:), 1)
        fz(jj) = sum(fz1(jj,:), 1)
      end do
      call vecA(xx, yy, zz, px, py, pz, fx, fy, fz, alphay, alphaz, aetax1, aetay1, aetaz1, &
                        aetax2, aetay2, aetaz2, nbath, nparticles, Axx, Ayy, Azz, Apx, Apy, Apz)
      call vecB_edges(dst, nparticles, dOmx, dOmy, dOmz)
      !call vecB_cool(dst, nparticles, dOmxc, dOmyc)
      xxi  = xx + Axx*dt
      yyi  = yy + Ayy*dt
      zzi  = zz + Azz*dt
      pxi = px + Apx*dt  + stermsBx*dOmx !+ stermsCx*dOmxc
      pyi = py + Apy*dt  + stermsBy*dOmy !+ stermsCy*dOmyc
      pzi = pz + Apz*dt  + stermsBz*dOmz !+ stermsCy*dOmyc
      fx = 0.0d0
      fy = 0.0d0
      fz = 0.0d0
      call coulombM(nparticles, xxi, yyi, zzi, fx2, fy2, fz2, invD2)
      do jj=1, nparticles, 1
        fx(jj) = sum(fx2(jj,:), 1)
        fy(jj) = sum(fy2(jj,:), 1)
        fz(jj) = sum(fz2(jj,:), 1)
      end do
      call vecA(xxi, yyi, zzi, pxi, pyi, pzi, fx, fy, fz, alphay, alphaz, aetax1, aetay1, aetaz1, &
                      aetax2, aetay2, aetaz2, nbath, nparticles, Axxi, Ayyi, Azzi, Apxi, Apyi, Apzi)
      xxN = xx + 0.5d0*(Axx + Axxi)*dt
      yyN = yy + 0.5d0*(Ayy + Ayyi)*dt
      zzN = zz + 0.5d0*(Azz + Azzi)*dt
      pxN = px + 0.5d0*(Apx + Apxi)*dt + stermsBx*dOmx !+ stermsCx*dOmxc
      pyN = py + 0.5d0*(Apy + Apyi)*dt + stermsBy*dOmy !+ stermsCy*dOmyc
      pzN = pz + 0.5d0*(Apz + Apzi)*dt + stermsBz*dOmz !+ stermsCy*dOmyc

      xx = xxN
      yy = yyN
      zz = zzN
      px = pxN
      py = pyN
      pz = pzN
      !print*,rank,  xx(1:5)

      tp = tp + dt

      if(tp .gt. tt) then
        !print*,"in we go"
        !print*, ps, save_freq
        if(mod(ps, save_freq) .eq. 0) then
          !print*, "saving"
          !xxS(:, saved) = xxS(:, saved) + xx
          !yyS(:, saved) = yyS(:, saved) + yy
          !zzS(:, saved) = zzS(:, saved) + zz
          pxS(:, saved) = pxS(:, saved) + px
          pyS(:, saved) = pyS(:, saved) + py
          pzS(:, saved) = pzS(:, saved) + pz
          px2S(:, saved) = px2S(:, saved) + px*px
          py2S(:, saved) = py2S(:, saved) + py*py
          pz2S(:, saved) = pz2S(:, saved) + pz*pz
          call local_energy(nparticles, alphay, alphaz, xx, yy, zz, invD1, px, py, pz, energy)
          call heat_current(nparticles, fx1, fy1, fz1, px, py, pz, hc)
          call current_Flux(hc, energy, xx, yy, zz, px, py, pz, nparticles, JJix, JJiy, JJiz)
          JJix_av = JJix_av + JJix
          JJiy_av = JJiy_av + JJiy
          JJiz_av = JJiz_av + JJiz
          !JJix_av_v(kk) = JJix_av_v(kk) + JJix/(nsteps-fin-1)
          !JJiy_av_v(kk) = JJiy_av_v(kk) + JJiy/(nsteps-fin-1)
          !JJiz_av_v(kk) = JJiz_av_v(kk) + JJiz/(nsteps-fin-1)
          saved = saved + 1
        end if
        if(saved .gt. nsave) continue = .False.
        !print*, continue
        ps = ps + 1
        if((mod(kk,10) .eq. 0)) then
          !px2S = px2S/kk
          !py2S = py2S/kk
          !pz2S = pz2S/kk
          !call mpi_reduce(xxS/kk, xxS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          !call mpi_reduce(yyS/kk, yyS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          !call mpi_reduce(zzS/kk, zzS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(px2S/kk, px2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(py2S/kk, py2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(pz2S/kk, pz2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(JJix/nsave, JJix_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(JJiy/nsave, JJiy_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(JJiz/nsave, JJiz_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          call mpi_reduce(kk, ntrajC, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
          !xxS_av  = xxS_av/procs
          !yyS_av  = yyS_av/procs
          !zzS_av  = zzS_av/procs
          px2S_av = px2S_av/procs
          py2S_av = py2S_av/procs
          pz2S_av = pz2S_av/procs
          JJix_av = JJix_av/ntrajC
          JJiy_av = JJiy_av/ntrajC
          JJiz_av = JJiz_av/ntrajC
          !print*, px2s_av(1:3,0)
          if(rank .eq. 0) then
            !open(unit=11, file="pos.dat")
            if(kk .eq. 1) then
              if(opened .eqv. .False.) then
                open(unit=11, file="pos.dat")
                opened = .True.
              end if
              !print*, yy(1:4)
              write(11,*)  xx(:), yy(:), zz(:)
            end if
            open(unit=12, file="psquared.dat")
            open(unit=13, file="heatflux.dat")
            !print*, xxS_av(1:4,:)
            do jj=1, nsave, 1
              !write(11,*) xxS_av(:,jj), yyS_av(:,jj), zzS_av(:,jj)
              write(12,*) px2S_av(:,jj), py2S_av(:,jj), pz2S_av(:,jj)
            end do
            write(13,*) JJix_av, JJiy_av, JJiz_av
            !close(unit=11)
            close(unit=12)
            close(unit=13)
          end if
        end if
      end if
    end do

    !JJix_av = JJix_av/nsave
    !JJiy_av = JJiy_av/nsave
    !JJiz_av = JJiz_av/nsave



    if(kk .eq. local_traj) then
      call mpi_barrier(mpi_comm_world, ierr)
      ! Last round of saving
      !px2S = pxS*pxS/(kk*kk)
      !py2S = pyS*pyS/(kk*kk)
      !pz2S = pzS*pzS/(kk*kk)
      !call mpi_reduce(xxS/kk, xxS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      !call mpi_reduce(yyS/kk, yyS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      !call mpi_reduce(zzS/kk, zzS_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(px2S/kk, px2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(py2S/kk, py2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(pz2S/kk, pz2S_av, nparticles*nsave, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(JJix/nsave, JJix_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(JJiy/nsave, JJiy_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(JJiz/nsave, JJiz_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(kk, ntrajC, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      xxS_av  = xxS_av/procs
      yyS_av  = yyS_av/procs
      zzS_av  = zzS_av/procs
      px2S_av = px2S_av/procs
      py2S_av = py2S_av/procs
      pz2S_av = pz2S_av/procs
      JJix_av = JJix_av/ntrajC
      JJiy_av = JJiy_av/ntrajC
      JJiz_av = JJiz_av/ntrajC
      if(rank .eq. 0) then
        !open(unit=11, file="pos.dat")
        open(unit=12, file="psquared.dat")
        open(unit=13, file="heatflux.dat")
        do jj=1, nsave, 1
          !write(11,*) xxS_av(:,jj), yyS_av(:,jj), zzS_av(:,jj)
          write(12,*) px2S_av(:,jj), py2S_av(:,jj), pz2S_av(:,jj)
        end do
        write(13,*) JJix_av, JJiy_av, JJiz_av
        close(unit=11)
        close(unit=12)
        close(unit=13)
      end if
    end if

  end do

  call mpi_finalize(ierr)



end program chainSolver

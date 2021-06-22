subroutine sp_seeds1(n_new_seeds, new_centers,  new_t0s, &
                     n_old_seeds, old_centers,  old_t0s, &
					 G, r_min, n_seeds, centers, t0s)
  ! This subroutine checks if seeds are born, or killed by other
  ! previously born particles. n_old seeds are already present in the computational
  ! volume.
  !
  ! n_new_seeds
  ! new_centers - array of dimension (n_new_seeds, 3)
  ! new_t0s - birth time, must be ordered! Array of dim (n_new_seeds)
  ! n_old_seeds
  ! old_centers - array of dimension (n_old_seeds, 3)
  ! old_t0s - birth time, must be ordered! Array of dim (n_old_seeds)
  ! G - growth speed. Scalar.
  ! r_min - distance from the edge of the spherulite, where already no new 
  !         spherulites are born
  !         This parameter helps to decrease the number of very small spherulites.
  !
  ! n_seeds
  ! centers (n_new_seeds + n_old_seeds, 3)
  ! t0s (n_new_seeds + n_old_seeds)
  !
  ! returns positions and birthtime of spontaneously born spheroulites
  !
  
  integer, intent(in) :: n_new_seeds, n_old_seeds
  real, dimension(n_new_seeds), intent(in) :: new_t0s
  real, dimension(n_new_seeds,3), intent(in) :: new_centers
  real, dimension(n_old_seeds), intent(in) :: old_t0s
  real, dimension(n_old_seeds,3), intent(in) :: old_centers

  real, intent(in) :: G
  real, intent(in) :: r_min
  
  integer, intent(out) :: n_seeds
  real, dimension(n_new_seeds+n_old_seeds), intent(out) :: t0s
  real, dimension(n_new_seeds+n_old_seeds,3), intent(out) :: centers

  integer :: it
  logical :: dead
  real, dimension(n_new_seeds+n_old_seeds) :: dd
  
  n_seeds = n_old_seeds
  centers(1:n_old_seeds,:) = old_centers(:,:)
  t0s(1:n_old_seeds) = old_t0s(:)

  do it = 1, n_new_seeds
    dead = .false.
    call arr_distance(new_centers(it,:), centers, n_seeds, dd(1:n_seeds))
    if (any(dd(1:n_seeds) < r_min + G*(new_t0s(it)-t0s))) then
      dead = .true.
    end if
     
!     do ip = 1, new_n_seeds
!        if (distance(new_centers(ip,:), centers(it,:)) < G*(t0s(it)-new_t0s(ip))) then
!           dead = .true.
!           exit
!        end if
!     end do
    
     if (.not. dead) then
        n_seeds = n_seeds + 1
        centers(n_seeds,:) = new_centers(it,:)
        t0s(n_seeds) = new_t0s(it)
     end if
  end do
end subroutine sp_seeds1

subroutine sp_seeds(n_new_seeds, new_centers,  new_t0s, &
			        G, r_min, n_seeds, centers, t0s)
  ! This subroutine checks if seeds are born, or killed by other
  ! previously born particles.
  !
  ! n_new_seeds
  ! new_centers - array of dimension (n_new_seeds, 3)
  ! new_t0s - birth time, must be ordered! Array of dim (n_new_seeds)
  ! G - growth speed. Scalar.
  ! r_min - distance from the edge of the spherulite, where already no new 
  !         spherulites are born
  !         This parameter helps to decrease the number of very small spherulites.
  !
  ! n_seeds
  ! centers (n_new_seeds + n_old_seeds, 3)
  ! t0s (n_new_seeds + n_old_seeds)
  !
  ! returns positions and birthtime of spontaneously born spheroulites
  !
  
  integer, intent(in) :: n_new_seeds
  real, dimension(n_new_seeds), intent(in) :: new_t0s
  real, dimension(n_new_seeds,3), intent(in) :: new_centers
  
  real, intent(in) :: G
  real, intent(in) :: r_min
  
  integer, intent(out) :: n_seeds
  real, dimension(n_new_seeds), intent(out) :: t0s
  real, dimension(n_new_seeds,3), intent(out) :: centers

  integer :: it
  logical :: dead
  real, dimension(n_new_seeds) :: dd
  
  n_seeds = 1
  centers(1,:) = new_centers(1,:)
  t0s(1) = new_t0s(1)

  do it = 2, n_new_seeds
    dead = .false.
    call arr_distance(new_centers(it,:), centers, n_seeds, dd(1:n_seeds))
    if (any(dd(1:n_seeds) < r_min + G*(new_t0s(it)-t0s))) then
      dead = .true.
    end if
     
!     do ip = 1, new_n_seeds
!        if (distance(new_centers(ip,:), centers(it,:)) < G*(t0s(it)-new_t0s(ip))) then
!           dead = .true.
!           exit
!        end if
!     end do
    
     if (.not. dead) then
        n_seeds = n_seeds + 1
        centers(n_seeds,:) = new_centers(it,:)
        t0s(n_seeds) = new_t0s(it)
     end if
  end do
end subroutine sp_seeds

subroutine sr_cdf(n_seeds, centers, t0s, LL, G, n_rnd_pos, t_conv, vols)
  ! calculates the cdf of a single realization, given the positions and 
  ! starting time of the seeds. 
  ! The subroutine assumes that all the seeds are born. Thus please run 
  ! sp_seeds subroutine before this one, if calculating spontaneous crystal growth.
  !
  ! n_seeds: number of spherulite seeds to grow
  ! centers: coordinates of the centers of spherulites (n_seeds,3)
  ! t0s: birthtime of the spherulites (n_seeds)
  ! LL: length of computational cube
  ! G: growing speed, scalar
  ! n_rnd_pos: number of random positions to check for crystallization time
  ! t_conv: to be returned, time of conversion at the positions (n_rnd_pos)
  ! vols: approximate volumes of the spherulites
  
  integer, intent(in) :: n_seeds, n_rnd_pos
  real, dimension(n_seeds, 3), intent(in) :: centers
  real, dimension(n_seeds), intent(in) :: t0s
  real, intent(in) :: LL, G
  real, dimension(n_rnd_pos), intent(out) :: t_conv
  real, dimension(n_seeds), intent(out) :: vols
  
  real, dimension(n_rnd_pos, 3) :: test_points
  real, dimension(n_seeds) :: t_conv_arr
  
  integer :: i, mloc

  real, dimension(n_seeds) :: dd_arr

  call init_random_seed()
  call random_number(test_points)
  test_points = test_points*LL
  vols = 0.0
  
  do i = 1, n_rnd_pos
    call arr_distance(test_points(i,:), centers, n_seeds, dd_arr)
    t_conv_arr = dd_arr/G+t0s
    mloc = minloc(t_conv_arr,1)
    t_conv(i) = t_conv_arr(mloc)
    vols(mloc) = vols(mloc)+1.0
    
    ! dd = distance(test_points(i,:),centers(1,:))
    ! t_conv_min = dd/G+t0s(1)
    ! do j = 2, n_seeds
    !    dd = distance(test_points(i,:),centers(j,:))
    !    tt = dd/G+t0s(j)
    !    if (tt < t_conv_min) then
    !       t_conv_min = tt
    !    end if
    ! end do
    ! t_conv(i) = t_conv_min

 end do

 vols = vols/n_rnd_pos*LL**3
end subroutine sr_cdf


function distance(c1,c2)
  real :: distance
  real, dimension(3) :: c1, c2
  
  distance = sqrt(sum((c1-c2)**2))
  return
end function distance


subroutine arr_distance(sc, arr, arr_length, dd)

  integer, intent(in) :: arr_length
  real, dimension(arr_length, 3), intent(in) :: arr
  real, dimension(3), intent(in) :: sc
  real, dimension(arr_length), intent(out) :: dd
  
  dd = sqrt(sum((arr-spread(sc,1,arr_length))**2, 2))
  
end subroutine arr_distance


subroutine init_random_seed()

  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
end subroutine init_random_seed


! subroutine mc_cdf_sp(n_seeds, LL, n_rnd_pos, G, t0s, t_min_arr)
!  !
!   ! calculates cummulative distribution function with spontaneous spherulite birth
!   ! t0s must be ordered!
!   ! 

!   integer, intent(in) :: n_seeds, n_rnd_pos
!   real, intent(in) :: LL, G
!   real, intent(out), dimension(n_rnd_pos) :: t_min_arr
!   real, intent(in), dimension(n_seeds) :: t0s

!   real, dimension(n_seeds, 3) :: centers, new_centers
!   real, dimension(n_seeds) :: new_t0s

!   integer :: ir, new_n_seeds

!   call init_random_seed()

!   call random_number(centers)
!   centers = centers*LL
     
!   ! call sp_seeds(n_seeds, centers, t0s, G, new_n_seeds, new_centers, new_t0s)
!   call sr_cdf(n_seeds, centers, new_t0s, n_rnd_pos, LL, G, t_min_arr)
  
! end subroutine mc_cdf_sp

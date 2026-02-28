


MODULE parameters
  
  ! atom in pixel data structure
  TYPE pixel_atom
    INTEGER                                 :: atom_id            ! the atom ID
    REAL, DIMENSION(3)                      :: coord              ! the atom coordinates on x, y and z
  END TYPE pixel_atom
  
  ! pixel data structure
  TYPE pixel
    INTEGER                                 :: pid                ! the pixel number
    INTEGER, DIMENSION(3)                   :: p_xyz              ! the pixel coordinates in the grid
    LOGICAL                                 :: tested             ! was the pixel checked already
    INTEGER                                 :: patoms             ! number of atom(s) in pixel
    TYPE(pixel_atom), DIMENSION(:), ALLOCATABLE   :: pix_atoms    ! list of atom(s) in pixel, to be allocated
    INTEGER                                 :: neighbors          ! number of neighbors for pixel
    INTEGER, DIMENSION(27)                  :: pixel_neighbors    ! the list of neighbor pixels, maximum 27
  END TYPE pixel
  
  ! pixel grid data structure
  TYPE pixel_grid
    INTEGER                                 :: pixels             ! total number of pixels in the grid
    INTEGER, DIMENSION(3)                   :: n_pix              ! number of pixel(s) on each axis
    INTEGER                                 :: n_xy               ! number of pixels in the plan xy
    TYPE (pixel), DIMENSION(:), ALLOCATABLE, TARGET :: pixel_list ! pointer to the pixels, to be allocated
  END TYPE pixel_grid
  
  ! distance data structure
  TYPE distance
    REAL                                    :: length             ! the distance in Angstrom squared
    REAL, DIMENSION(3)                      :: Rij                ! vector components of x, y and z
  END TYPE distance
  
  !
  ! the following are considered to be provided by the user
  !

  ! model description
  INTEGER                                   :: atoms              ! the total number of atom(s)
  REAL, DIMENSION(atoms,3)                  :: c_coord            ! list of Cartesian coordinates
  REAL                                      :: cutoff             ! the cutoff to define first neighbors
  REAL                                      :: cutoff_squared     ! squared value of the cutoff to define first neighbors
  
  ! model box description
  REAL, DIMENSION(3)                        :: l_params           ! lattice a, b and c
  REAL, DIMENSION(3,3)                      :: cart_to_frac       ! Cartesian to fractional coordinates matrix
  REAL, DIMENSION(3,3)                      :: frac_to_cart       ! fractional to Cartesian coordinates matrix

END MODULE parameters


SUBROUTINE set_pbc_shift (grid, pixel_coord, pbc_shift)
  
  USE parameters
  IMPLICIT NONE
  
  TYPE (pixel_grid), INTENT(IN)            :: grid                ! the pixel grid
  INTEGER, DIMENSION(3), INTENT(IN)        :: pixel_coord         ! the pixel coordinates in the grid
  INTEGER, DIMENSION(3,3,3), INTENT(INOUT) :: pbc_shift           ! the shift, correction, to be calculated
  
  pbc_shift(:,:,:) = 0                                            ! initialization without any shift

  if ( pixel_coord(1) .eq. 1 ) then                               ! pixel position on 'x' is min
    pbc_shift(1,:,:) = grid%n_pix(1)
  else if ( pixel_coord(1) .eq. grid%n_pix(1) ) then              ! pixel position on 'x' is max
    pbc_shift(3,:,:) = - grid%n_pix(1)
  endif
  
  if ( pixel_coord(2) .eq. 1 ) then                               ! pixel position on 'y' is min
    pbc_shift(:,1,:) = pbc_shift(:,1,:) + grid%n_xy
  else if ( pixel_coord(2) .eq. grid%n_pix(2) ) then              ! pixel position on 'y' is max
    pbc_shift(:,3,:) = pbc_shift(:,3,:) - grid%n_xy
  endif
  
  if ( pixel_coord(3) .eq. 1 ) then                               ! pixel position on 'z' is min
    pbc_shift(:,:,1) = pbc_shift(:,:,1) + grid%pixels
  else if ( pixel_coord(3) .eq. grid%n_pix(3) ) then              ! pixel position on 'z' is max
    pbc_shift(:,:,3) = pbc_shift(:,:,3) - grid%pixels
  endif
  
END SUBROUTINE set_pbc_shift


SUBROUTINE add_atom_to_pixel (the_pixel, pixel_coord, atom_id, atom_coord)

  USE parameters
  IMPLICIT NONE

  TYPE (pixel), INTENT(INOUT) :: the_pixel          ! pointer to the pixel
  INTEGER, DIMENSION(3), INTENT(IN) :: pixel_coord  ! the pixel coordinates in the grid
  INTEGER, INTENT(IN)                :: atom_id     ! the atom ID number
  INTEGER                            :: axis        ! loop iterator axis id (1=x, 2=y, 3=z)
  REAL, DIMENSION(3), INTENT(IN)     :: atom_coord  ! the atomic coordinates

  if (the_pixel%patoms .eq. 0) then
    ! if the pixel do not contains any atom yet, then save its coordinates in the grid
    do axis = 1 , 3
      the_pixel%p_xyz(axis) = pixel_coord(axis)
    enddo
     ! allocate the memory to store the first pixel_atom information
    allocate(the_pixel%pix_atoms(1))
  else
    ! otherwise reallocate memory to store the new pixel_atom informatio
    the_pixel%pix_atoms = realloc(the_pixel%pix_atoms, the_pixel%patoms+1)
  endif
  ! increment the number of atom(s) in the pixel
  the_pixel%patoms = the_pixel%patoms + 1
  the_pixel%pix_atoms(the_pixel%patoms)%atom_id = atom_id
  do axis = 1 , 3
    the_pixel%pix_atoms(the_pixel%patoms)%coord(axis) = atom_coord(axis)
  enddo

END SUBROUTINE



! Preparation of the pixel grid
SUBROUTINE prepare_pixel_grid (use_pbc, grid)
  
  USE parameters
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN)              :: use_pbc          ! flag to set if PBC are used or not
  TYPE (pixel_grid), INTENT(INOUT) :: grid             ! the pixel grid to prepare
  INTEGER                          :: axis             ! loop iterator axis id (1=x, 2=y, 3=z)
  INTEGER                          :: aid              ! loop iterator atom number (1, atoms)
  INTEGER                          :: pixel_num        ! pixel number in the grid
  INTEGER, DIMENSION(3)            :: pixel_pos        ! pixel coordinates in the grid
  REAL, DIMENSION(3)               :: cmin, cmax       ! real precision coordinates min, max
  REAL, DIMENSION(3)               :: f_coord          ! real precision fractional coordinates
  
  if ( .not. use_pbc ) then                            ! without periodic boundary conditions
    do axis = 1 , 3
      cmin(axis) = c_coord(1,axis)
      cmax(axis) = c_coord(1,axis)
    enddo
    do aid = 2 , atoms                                 ! for all atoms
      do axis = 1 , 3                                  ! for x, y and z
        cmin(axis) = min(cmin(axis), c_coord(aid,axis))
        cmax(axis) = max(cmax(axis), c_coord(aid,axis))
      enddo
    enddo
    do axis = 1 , 3                                    ! for x, y and z
      ! Number of pixel(s) on axis 'axis'
      grid%n_pix(axis) = INT((cmax(axis) - cmin(axis)) / cutoff) + 1
    enddo
  else                                                 ! using periodic boundary conditions
    do axis = 1 , 3                                    ! for x, y and z
      ! Number of pixel(s) on axis 'axis'
      grid%n_pix(axis) = INT(l_params(axis) / cutoff) + 1
    enddo
  endif
  do axis = 1 , 3                                      ! for x, y and z
    ! correction if the number of pixel(s) on 'axis' is too small
    if ( grid%n_pix(axis) .lt. 4 ) then
      grid%n_pix(axis) = 1
    endif
  enddo
  
  grid%n_xy = grid%n_pix(1) * grid%n_pix(2)            ! Number of pixels on the plan 'xy'
  grid%pixels = grid%n_xy * grid%n_pix(3)              ! Total number of pixels in the grid

  allocate (grid%pixel_list(grid%pixels))
  do pixel_num = 1 , grid%pixels
    grid%pixel_list(pixel_num)%pid = pixel_num
    grid%pixel_list(pixel_num)%patoms = 0
    grid%pixel_list(pixel_num)%tested = .false.
  enddo
  if ( .not. use_pbc ) then                            ! without periodic boundary conditions
    do aid = 1 , atoms                                 ! for all atoms
      do axis = 1 , 3                                  ! for x, y and z
        pixel_pos(axis) = INT( (c_coord(aid,axis) - cmin(axis) )/ cutoff)
      enddo
      pixel_num = pixel_pos(1) + pixel_pos(2) * grid%n_pix(1) + pixel_pos(3) * grid%n_xy + 1
      call add_atom_to_pixel (grid%pixel_list(pixel_num), pixel_pos, aid, c_coord(aid,:))
    enddo
  else                                                 ! using periodic boundary conditions
    do aid = 1 , atoms                                 ! for all atoms
     f_coord = MATMUL ( c_coord(aid), cart_to_frac )
     do axis = 1 , 3                                   ! for x, y and z
       f_coord(axis) = f_coord(axis) - floor(f_coord(axis))
       pixel_pos(axis) = INT(f_coord(axis) * grid%n_pix(axis))
     enddo
     pixel_num = pixel_pos(1) + pixel_pos(2) * grid%n_pix(1) + pixel_pos(3) * grid%n_xy + 1
     call add_atom_to_pixel (grid%pixel_list(pixel_num), pixel_pos, aid, f_coord)
    enddo
  endif

END SUBROUTINE prepare_pixel_grid


SUBROUTINE find_pixel_neighbors (use_pbc, the_grid, the_pix)
  
  USE parameters
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN)                  :: use_pbc                 ! flag to set if PBC are used or not
  TYPE (pixel_grid), INTENT(INOUT)     :: the_grid                ! pointer to the pixel grid
  TYPE (pixel), POINTER, INTENT(INOUT) :: the_pix                 ! pointer to the pixel
  INTEGER                              :: axis                    ! loop iterator axis id (1=x, 2=y, 3=z)
  INTEGER                              :: xpos, ypos, zpos        ! neighbor position on x, y and z
  INTEGER, DIMENSION(3)                :: l_start = (\1, 1, 1\)   ! loop iterators starting value
  INTEGER, DIMENSION(3)                :: l_end = (\3, 3, 3\)     ! loop iterators ending value
  INTEGER, DIMENSION(3)                :: pmod = (\-1, 0, 1\)     ! position modifiers
  INTEGER                              :: nnp                     ! number of neighbors for pixel
  INTEGER                              :: nid                     ! neighbor id for pixel
  INTEGER, DIMENSION(3,3,3)            :: pbc_shift               ! shift, correction, due to PBC
  LOGICAL                              :: boundary=.false.        ! is pixel on the boundary of the grid
  LOGICAL                              :: keep_neighbor = .true.  ! keep or not neighbor during analysis
  
  if ( use_pbc ) then
    call set_pbc_shift (the_grid, the_pix%p_xyz, pbc_shift)
  else
    do axis = 1 , 3
      if ( the_pix%p_xyz(axis) .eq. 1 .or. the_pix%p_xyz(axis) .eq. the_grid%n_pix(axis) ) then
        boundary = .true.
      endif
    enddo
  endif
  do axis = 1 , 3
    if ( the_grid%n_pix(axis) .eq. 1 ) then
      l_start(axis) = 2
      l_end(axis) = 2
    endif
  enddo
  nnp = 0
  do xpos = l_start(1) , l_end(1)
    do ypos = l_start(2) , l_end(2)
      do zpos = l_start(3) , l_end(3)
        keep_neighbor = .true.
        if ( .not. use_pbc .and. boundary ) then
          if ( the_pix%p_xyz(1) .eq. 1 .and. xpos .eq. 1 ) then
            keep_neighbor = .false.
          else if ( the_pix%p_xyz(1) .eq. the_grid%n_pix(1) .and. xpos .eq. 3 ) then
            keep_neighbor = .false.
          else if ( the_pix%p_xyz(2) .eq. 1 .and. ypos .eq. 1 ) then
            keep_neighbor = .false.
          else if ( the_pix%p_xyz(2) .eq. the_grid%n_pix(2) .and. ypos .eq. 3 ) then
            keep_neighbor = .false.
          else if ( the_pix%p_xyz(3) .eq. 1 .and. zpos .eq. 1 ) then
            keep_neighbor = .false.
          else if ( the_pix%p_xyz(3) .eq. the_grid%n_pix(3) .and. zpos .eq. 3 ) then
            keep_neighbor = .false.
          endif
        endif
        if ( keep_neighbor ) then
          ! evaluating neighbor pixel number in the grid
          nid = the_pix%pid + pmod(xpos) + pmod(ypos) * the_grid%n_pix(1) + pmod(zpos) * the_grid%n_xy
          if ( use_pbc ) then
            ! correcting the value if PBC are used
            nid = nid + pbc_shift(xpos, ypos, zpos)
          endif
          the_pix%pixel_neighbors(nnp) = nid
          nnp = nnp + 1
        endif
      enddo
    enddo
  enddo
  the_pix%neighbors = nnp
  
END SUBROUTINE find_pixel_neighbors



! evaluating the interatomic distance between 2 pixel atoms
SUBROUTINE evaluate_distance (use_pbc, at_i, at_j, dist)
  
  USE parameters
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN)                    :: use_pbc   ! flag to set if PBC are used or not
  TYPE (pixel_atom), POINTER, INTENT(IN) :: at_i      ! first pixel atom
  TYPE (pixel_atom), POINTER, INTENT(IN) :: at_j      ! second pixel atom
  TYPE (distance), INTENT(INOUT)         :: dist      ! calculation results
  INTEGER                                :: axis      ! loop iterator axis id (1=x, 2=y, 3=z)
  
  do axis = 1 , 3
    dist%Rij(axis) = at_i%coord(axis) - at_j%coord(axis)
  enddo
  if ( use_pbc ) then
    ! then the pixel_atom's coordinates are in corrected fractional format
    do axis = 1 , 3
      dist%Rij(axis) = dist%Rij(axis) - AnINT(dist%Rij(axis))
    enddo
    ! transform back to Cartesian coordinates
    dist%Rij = MATMUL( dist%Rij, frac_to_cart )
  endif
  
  dist%length = 0.0
  do axis = 1 , 3
    dist%length = dist%length + dist%Rij(axis) * dist%Rij(axis)
  enddo
  ! returning the 'distance' data structure that contains:
  ! - the squared value for Dij: no time consuming square root calculation !
  ! - the components of the distance vector on x, y and z
END SUBROUTINE evaluate_distance



SUBROUTINE pixel_search_for_neighbors (use_pbc)

  USE parameters
  IMPLICIT NONE

  LOGICAL, INTENT(IN)   :: use_pbc                      ! flag to set if PBC are used or not
  TYPE (pixel_grid)     :: all_pixels                   ! the pixel grid to analyze
  INTEGER               :: pix, pjx                     ! integer pixel ID numbers
  INTEGER               :: aid, bid                     ! integer loop atom numbers
  INTEGER               :: l_start, l_end               ! integer loop modifier
  INTEGER               :: pid
  TYPE (pixel_atom), POINTER :: pix_i, pix_j            ! pointers on pixel data structure
  TYPE (atom), POINTER  :: at_i, at_j                   ! pointers on pixel_atom data structure
  TYPE (distance)       :: Dij                          ! distance data structure
  
  call prepare_pixel_grid (use_pbc, all_pixels)
  ! note that the pixel grid 'all_pixels' must be prepared before the following
  ! for all pixels in the grid
  do pix = 1 , all_pixels%pixels
    ! setting 'pix_i' as pointer to pixel number 'pix'
    pix_i => all_pixels%pixel_list(pix)
    ! if pixel 'pix_i' contains atom(s)
    if ( pix_i%patoms .gt. 0 ) then

      ! search for pixel neighbbors of 'pix_'
      call find_pixel_neighbors (use_pbc, all_pixels, pix_i)

      ! testing all 'pix_i' neighbor pixels
      do pid = 1 , pix_i%neighbors
        pjx = pix_i%pixel_neighbors(pid)
        ! setting 'pix_j' as pointer to pixel number 'pjx'
        pix_j => all_pixels%pixel_list(pjx)
        ! checking pixel 'pix_j' if it:
        ! - contains atom(s)
        ! - was not tested, otherwise the analysis would have been performed already
        if ( pix_j%patoms .gt. 0 .and. .not. pix_j%tested ) then
          ! If 'pix_i' and 'pix_j' are the same, only test pair of different atoms
          if ( pjx .eq. pix ) then
            l_end = 1
          else
            l_end = 0
          endif
          ! for all atom(s) in 'pix_i'
          do aid = 1 , pix_i%patoms - l_end
            ! set pointers to the first atom to test
            at_i => pix_i%pix_atoms(aid)
            if ( pjx .eq. pix ) then
              l_start = aid + 1
            else
              l_start = 1
            endif
            ! for all atom(s) in 'pix_j'
            do bid = l_start , pix_j%patoms
              ! set pointers to the second atom to test
              at_j => pix_j%pix_atoms(bid)
              ! evaluate interatomic distance
              call evaluate_distance (use_pbc, at_i, at_j, Dij)
              if ( Dij%length .lt. cutoff_squared ) then
                ! this is a first neighbor bond !
              endif
            enddo
          enddo
        endif
      enddo
      ! store that pixel 'pix_i' was tested
      pix_i%tested = .true.
    endif
  enddo
  
END SUBROUTINE pixel_search_for_neighbors

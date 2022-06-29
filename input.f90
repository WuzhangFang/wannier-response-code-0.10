module input
use constants
use iso_fortran_env, only: input_unit
implicit none
!--------------------------------------------------------------------------
!Contains all the necessary input parameters and subroutine for reading them.
!--------------------------------------------------------------------------

type input_type
    character(len=100) :: linres
    character(len=100) :: model
    integer :: nk(3)
    real(dp), allocatable :: gam(:)
    real(dp) :: Ef
    integer :: n_even
    integer, allocatable :: for_even(:)
    integer :: n_odd
    integer, allocatable :: for_odd(:)
    real(dp) :: eps
    integer :: n_projs
    integer, allocatable :: atoms(:)
    integer :: n_op1
end type

type projs_type
    integer :: n_atoms
    integer, allocatable :: n_orbs(:)
    integer, allocatable :: orbs(:,:)
end type

type bands_input_type
    integer :: n_kpts
    integer :: n_sym_kpts
    real(dp), allocatable :: sym_kpts(:,:)
    character(1), allocatable :: kpts_names(:)
    logical :: print_Ham
end type

contains

subroutine read_input(inp)
    type(input_type), intent(out) :: inp

    integer :: ngam
    real(dp) :: gam1, gamf
    integer :: i
    integer :: u

    open(newunit=u, file=input_file, status='old')

    read(u,*) inp%linres
    read(u,*) inp%model
    read(u,*) inp%nk
    read(u,*) ngam,gam1,gamf
    read(u,*) inp%Ef
    read(u,*) inp%n_even
    allocate(inp%for_even(inp%n_even))
    if (inp%n_even /= 0) read(u,*) inp%for_even
    read(u,*) inp%n_odd
    allocate(inp%for_odd(inp%n_odd))
    if (inp%n_odd /= 0) read(u,*) inp%for_odd
    read(u,*) inp%eps
    read(u,*) inp%n_projs

    if ( inp%n_projs > 0 .and. ( inp%linres == 'cond' .or. inp%linres == 'she') ) then
        write(*,*) 'Projections only allowed for cisp.'
        write(*,*) 'setting n_projs = 0'
        inp%n_projs = 0
    end if

    if (inp%n_projs > 0 ) then
        allocate(inp%atoms(inp%n_projs))
        read(u,*) inp%atoms
    end if
    
    allocate(inp%gam(ngam))
    if (ngam == 1) then
        inp%gam(1) = gam1
    else
        do i=1,ngam
            inp%gam(i) = 10._dp ** ( (log10(gamf)-log10(gam1)) * (i-1)/(ngam-1) + log10(gam1))
        end do
    endif

    if (inp%linres == 'cisp') then
        inp%n_op1 = max(1,inp%n_projs)*3
    else if (inp%linres == 'cond') then
        inp%n_op1 = 3
    else if (inp%linres == 'she') then
        inp%n_op1 = 9
    else
        write(*,*) inp%linres, 'is a wrong switch'
        stop 'wrong switch'
    endif 

end subroutine

subroutine read_projs(projs)
    type(projs_type), intent(out) :: projs

    integer :: i
    integer :: max_norb
    integer :: u

    open(newunit=u, file=projs_file, status='old')

    read(u,*) projs%n_atoms 
    allocate(projs%n_orbs(projs%n_atoms))
    do i=1,projs%n_atoms
        read(u,*) projs%n_orbs(i)
        read(u,*) 
    end do
    max_norb = maxval(projs%n_orbs)
    allocate(projs%orbs(max_norb,projs%n_atoms))
    rewind(u)
    read(u,*)
    do i=1,projs%n_atoms
        read(u,*) 
        read(u,*) projs%orbs(:,i)
    end do

end subroutine

subroutine read_bands_input(band_inp)

    type(bands_input_type), intent(out) :: band_inp

    integer :: i,u,iostatus

    open(newunit=u, file=bands_input_file, status='old')

    read(u,*) band_inp%n_kpts
    read(u,*) band_inp%n_sym_kpts

    allocate(band_inp%sym_kpts(3,band_inp%n_sym_kpts))
    allocate(band_inp%kpts_names(band_inp%n_sym_kpts))

    do i=1,band_inp%n_sym_kpts
        read(u,*) band_inp%kpts_names(i), band_inp%sym_kpts(:,i)
    end do

    read(u,*,IOSTAT=iostatus) band_inp%print_Ham
    if ( iostatus /= 0 ) then
        band_inp%print_Ham = .false.
    end if

end subroutine

subroutine read_structure(latt_vecs)
    !--------------------------------------------------------------------------
    !Reads the lattice vectors from a POSCAR file
    !--------------------------------------------------------------------------

    real(dp) :: latt_vecs(3,3)

    real(dp) :: latt_const
    integer :: u

    open(newunit=u, file=struct_file, status='old')

    !skips the first line
    read(u,*)
    read(u,*) latt_const
    !reads the lattice vectors
    read(u,*) latt_vecs(:,1)
    read(u,*) latt_vecs(:,2)
    read(u,*) latt_vecs(:,3)

    latt_vecs = latt_vecs * latt_const
    
    close(u)

end subroutine

end module




module wann_model_mod
!--------------------------------------------------------------------------
!This is a module for working with tight-binding Hamiltonians produced by
!Wannier90.
!It reads the wannier90_hr.dat file using the init routine and contains 
!the subroutines necessary for linear response calculations:
!   create_Hk: creates the k-dependent Hamiltonian at point k.
!   create_vk: constructs the velocity operator at point k
!   create_smat: creates the spin-operator
!--------------------------------------------------------------------------

use generic_model
use constants
use input
use blas95

type, extends (gen_model) :: wann_model

    integer :: n_r
    integer, allocatable :: r_deg(:)
    integer, allocatable :: rpts(:,:)
    complex(dp), allocatable :: Hr(:,:,:)

    contains
        procedure :: init
        procedure :: create_smat
        procedure :: create_Hk
        procedure :: create_vk
end type

contains


subroutine init(self)
    !--------------------------------------------------------------------------
    !Reads the hamiltonian matrix elements in the Wannier basis outputted by
    !wannier 90.
    !
    !Returns  (module variables):
    !   n_orb (integer): number of Wannier functions used
    !   n_r (integer): number of Wigner-Seitz grid points
    !   r_deg (integer(:)): degeneracies of the Wigner-Seitz grid points
    !   rpts (integer(:,:)): Wigner-Seitz grid points
    !   Hr (complex(:,:,:)): the Hamiltonian matrix elements stored as:
    !       the first index is the point number (as in rpts), the other two are
    !       orbital numbers
    !--------------------------------------------------------------------------

    class(wann_model) :: self

    integer :: u
    integer :: i
    integer :: nl
    integer :: stat
    integer :: n_wann
    integer :: r(3),r_curr(3)
    integer :: r_i
    integer :: orb1,orb2
    real(dp) :: val(2)

    open(newunit=u, file=Hr_file, status='old')

    read(u,*) 
    read(u,*) self%n_orb
    read(u,*) self%n_r

    n_wann = self%n_orb

    allocate(self%r_deg(self%n_r))
    allocate(self%rpts(3,self%n_r))
    allocate(self%Hr(self%n_r,n_wann,n_wann))

    self%Hr = 0._dp

    !read the degeneracies of the Wigner-Seitz grid points
    nl = self%n_r/15
    do i=1,nl
        read(u,*) self%r_deg(1+(i-1)*15:i*15)
    end do
    read(u,*) self%r_deg(nl*15+1:self%n_r)

    r_i = 0
    !reads the actual Hamiltonian and makes a list of all positions
    do 
        read(u,*,iostat=stat) r(:), orb1, orb2, val(:)
        if (stat /= 0) exit

        if (r_i == 0) then 
            r_i = 1
            r_curr = r
            self%rpts(:,r_i) = r
        endif 

        if (.not. all(r == r_curr)) then
            r_i = r_i + 1
            self%rpts(:,r_i) = r
            r_curr = r
        endif

        self%Hr(r_i,orb1,orb2) = cmplx(val(1),val(2))

    end do

    close(u)

end subroutine

subroutine create_Hk(self,k,Hk)
    !--------------------------------------------------------------------------
    !Creates the k-dependent Hamiltonian at point k by Fourrier transforming Hr
    !
    !Description:
    !   H_k(n,m) = sum_R exp(i*k*R)H_R(n,m)
    !
    !Notes:
    !   !!!!Only the lower diagonal part of Hk is calculated!!!!
    !
    !Args:
    !   r_deg: degeneracies of the Weigner-Seitz grid points 
    !   rpts: Weigner-Seitz grid points
    !   Hr: the Hamiltonian in the basis of Wannier functions
    !   k: k point
    !
    !Returns:
    !   Hk: the lower diagonal part of Hk
    !--------------------------------------------------------------------------

    class(wann_model) :: self
    real(dp), intent(in) :: k(3)
    complex(dp), allocatable, intent(out) :: Hk(:,:)

    integer :: n_r !number of R points
    integer :: n_wann
    integer :: i,n,m
    complex(dp), allocatable :: exps_vec(:)
    real(dp) :: latt_vecs(3,3) !lattice vectors
    real(dp) :: r(3)

    n_r = size(self%Hr,1)
    n_wann = self%n_orb

    allocate(Hk(n_wann,n_wann))
    Hk(:,:) = 0_dp

    allocate(exps_vec(n_r))

    call read_structure(latt_vecs)

    !to speed up the calculation we first calculate exp(i*k*R) and then calculate the fourrier
    !transform as a dot product
    do i=1,self%n_r
        !express R in cartesian coordinates in the units of Ang
        r = self%rpts(1,i) * latt_vecs(:,1) + self%rpts(2,i) * latt_vecs(:,2) + self%rpts(3,i) * latt_vecs(:,3) 
        !the conjugate is here becase dotc also contains a conjg and we dont want that
        !has to be divided by r_deg, don't really know why
        exps_vec(i) = conjg(exp(ii*dot_product(r,k)) / self%r_deg(i))
    end do

    !only the lower diagonal part is calculated
    do m=1,n_wann
        do n=m, n_wann
            
            Hk(n,m) = dotc(exps_vec,self%Hr(:,n,m))

        end do
    end do

end subroutine

subroutine create_vk(self,k,vk,ind)
    !--------------------------------------------------------------------------
    !Creates the velocity operator at point k as a derivative of Hk
    !
    !Description:
    !   Very similar to how Hk is calculated
    !   v_k(n,m) = 1/hbar * \partial(Hk)\partial(k) = 
    !       = 1/hbar * sum_R exp(i*k*R)H_R(n,m) * i * R
    !   When ind is present only the ind component of the velocity operator
    !   is calculated otherwise all three components are calculated
    !
    !Args: (needs output from read_Hr)
    !   r_deg: Degeneracies of Weigner-Seitz grid points
    !   rpts: Weigner-Seitz grid points
    !   Hr: the Wannier Hamiltonian
    !   k: k point 
    !   ind (optional): when present, only this component of the v operator is 
    !       returned
    !
    !Returns:
    !   vk: the velocity operator at point k
    !--------------------------------------------------------------------------

    class(wann_model) :: self
    real(dp), intent(in) :: k(3)
    complex(dp), allocatable, intent(out) :: vk(:,:,:)
    integer, intent(in), optional :: ind

    integer :: n_r
    integer :: n_wann
    integer :: i,j,n,m
    complex(dp), allocatable :: exps_vec(:)
    real(dp) :: latt_vecs(3,3)
    real(dp) :: r(3)

    n_r = size(self%Hr,1)
    n_wann = self%n_orb

    allocate(vk(n_wann,n_wann,3))
    vk(:,:,:) = 0_dp
    allocate(exps_vec(n_r))

    call read_structure(latt_vecs)

    !create the upper diagonal part of vk
    do j=1,3
        if (present(ind)) then !if ind is present and j/=ind skip to the next cycle
            if (ind /= j) cycle
        endif
        do i = 1,n_r
            r = self%rpts(1,i) * latt_vecs(:,1) + self%rpts(2,i) * latt_vecs(:,2) + self%rpts(3,i) * latt_vecs(:,3)
            exps_vec(i) = conjg(exp(ii*dot_product(r,k)) * ii * r(j) / self%r_deg(i))
        end do
        do m=1,n_wann
            do n=m, n_wann
                vk(n,m,j) = dotc(exps_vec,self%Hr(:,n,m))
            end do
        end do

        !use the fact the vk must be Hermitian to find the lower part
        do m=2,n_wann
            do n=1,m-1
                vk(n,m,j) = conjg(vk(m,n,j)) 
            end do
        end do
    end do

    vk = vk / hbar

end subroutine

subroutine create_smat(self,smat,atom,projs)
    !--------------------------------------------------------------------------
    !Defines the matrix of the spin-operator
    !
    !Notes: 
    !   Assumes that first half othe basis functions are spin-up and the other half
    !   are spin-down.
    !   dimensionless, i.e. S=sigma. Divide by 2 to get spin in the units of hbar
    !
    !Args:
    !   n_wann: magnitude(dimension) of the Hamiltonian
    !   atom(optional): number of atom on which the spin-operator is projected
    !       if atom is not present, no projection is done
    !   projs(optional): must be present if atom is present. defines how the 
    !       wannier orbitals are localized, which defines how the projections are
    !       done. This is outputed by read_projs
    !
    !Returns:
    !   smat: all three components of the Hamiltonian
    !--------------------------------------------------------------------------
    class(wann_model) :: self
    complex(dp), allocatable, intent(out) :: smat(:,:,:)
    integer, optional, intent(in) :: atom
    type(projs_type), optional, intent(in) :: projs

    integer :: i,orb,n_wann

    n_wann = self%n_orb

    allocate(smat(n_wann,n_wann,3))
    smat(:,:,:) = 0._dp
    
    ! no projection
    if ((.not. present(atom)) .or. (atom == 0)) then
    
        !sigma_x
        do i =1,n_wann/2
            smat(i,n_wann/2+i,1) = 1._dp
            smat(n_wann/2+i,i,1) = 1._dp
        end do

        !sigma_y
        do i =1,n_wann/2
            smat(i,n_wann/2+i,2) = (0._dp,-1._dp)
            smat(n_wann/2+i,i,2) = (0._dp,1._dp)
        end do

        !sigma_z
        do i =1,n_wann/2
            smat(i,i,3) = 1._dp 
            smat(n_wann/2+i,n_wann/2+i,3) = -1._dp 
        end do
    ! projection on # atom
    else
        !sigma_x
        do i =1,projs%n_orbs(atom)/2
            orb = projs%orbs(i,atom)
            smat(orb,n_wann/2+orb,1) = 1._dp
            smat(n_wann/2+orb,orb,1) = 1._dp
        end do

        !sigma_y
        do i =1,projs%n_orbs(atom)/2
            orb = projs%orbs(i,atom)
            smat(orb,n_wann/2+orb,2) = (0._dp,-1._dp)
            smat(n_wann/2+orb,orb,2) = (0._dp,1._dp)
        end do

        !sigma_z
        do i =1,projs%n_orbs(atom)/2
            orb = projs%orbs(i,atom)
            smat(orb,orb,3) = 1._dp 
            smat(n_wann/2+orb,n_wann/2+orb,3) = -1._dp 
        end do
    endif

end subroutine

end module


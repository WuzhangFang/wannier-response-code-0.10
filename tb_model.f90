module tb_model_mod
!--------------------------------------------------------------------------
!This module reads a tight-binding Hamiltonian produced by python sd_model.py
!script and provides routines for linear response calculation. The routines
!are similar to the wann_model and have the same output but work a little 
!differently.
!--------------------------------------------------------------------------

use generic_model
use constants
use input
use blas95

type, extends (gen_model) :: tb_model

    integer :: n_Hrs
    real(dp), allocatable :: rhos(:,:)
    integer, allocatable :: orbs(:,:)
    complex(dp), allocatable :: Hr(:)

    contains
        procedure :: init
        procedure :: create_smat
        procedure :: create_Hk
        procedure :: create_vk
end type

contains

subroutine init(self)
    !--------------------------------------------------------------------------
    !Reads the hamiltonian matrix elements in the tight-binding basis outputted by
    !sd_model.py.
    !
    !Returns (module variables):
    !   n_orb (integer): number of basis functions
    !   n_Hrs (integer): number of Hamiltonian matrix elements
    !   rhos (real(3,n_orb)): rhos that appear in the exponential for each matrix element
    !       rho is given in fractional coordinates and is defined as:
    !   orbs (integer(3,n_orb): the orbitals of each matrix element
    !   Hr (complex(:)): the Hamiltonian matrix elements stored as:
    !       Hr(i) = < orbs(1,i)|Hr|orbs(2,i) >
    !--------------------------------------------------------------------------

    class(tb_model) :: self

    integer :: u
    integer :: i
    integer :: n_wann
    real(dp) :: r(3)
    integer :: orb1,orb2
    real(dp) :: val(2)

    open(newunit=u, file='tb_hr.dat', status='old')

    read(u,*) 
    read(u,*) self%n_orb
    read(u,*) self%n_Hrs

    n_wann = self%n_orb

    allocate(self%rhos(3,self%n_Hrs))
    allocate(self%orbs(2,self%n_Hrs))
    allocate(self%Hr(self%n_Hrs))

    self%Hr = 0._dp

    do i=1,self%n_Hrs
        read(u,*) r(:), orb1, orb2, val(:)
        self%rhos(:,i) = r(:)
        self%orbs(1,i) = orb1
        self%orbs(2,i) = orb2
        self%Hr(i) = cmplx(val(1),val(2))
    end do

    close(u)

end subroutine

subroutine create_Hk(self,k,Hk)
    !--------------------------------------------------------------------------
    !Creates the k-dependent Hamiltonian at point k by Fourrier transforming Hr
    !
    !Description:
    !   H_k(n,m) = sum_R exp(i*k*rho)H_R(n,m)
    !   rho = R + p_m - p_n,
    !   we consider that orbital n is in the first unit cell, R is then vector
    !   of the cell of orbital m
    !   p_n,p_m are coordinates of atoms on which orbitals n,m are defined in
    !   the first unit cell
    !
    !Args:
    !   k (real(3)): the k point
    !Returns:
    !   Hk: the lower diagonal part of Hk
    !--------------------------------------------------------------------------

    class(tb_model) :: self
    real(dp), intent(in) :: k(3)
    complex(dp), allocatable, intent(out) :: Hk(:,:)

    integer :: n_wann
    integer :: i,n,m
    real(dp) :: latt_vecs(3,3) !lattice vectors
    real(dp) :: r(3)

    n_wann = self%n_orb

    allocate(Hk(n_wann,n_wann))
    Hk(:,:) = 0_dp

    call read_structure(latt_vecs)

    do i=1,self%n_Hrs
        !express rho in cartesian coordinates in the units of Ang
        r = self%rhos(1,i) * latt_vecs(:,1) + self%rhos(2,i) * latt_vecs(:,2) + self%rhos(3,i) * latt_vecs(:,3) 
        n = self%orbs(1,i)
        m = self%orbs(2,i)
        Hk(n,m) = Hk(n,m) + exp(ii*dot_product(k,r)) * self%Hr(i)
    end do

end subroutine

subroutine create_vk(self,k,vk,ind)
    !--------------------------------------------------------------------------
    !Creates the velocity operator at point k as a derivative of Hk
    !
    !Description:
    !   Very similar to how Hk is calculated
    !   v_k(n,m) = 1/hbar * \partial(Hk)\partial(k) = 
    !       = 1/hbar * sum_R exp(i*k*rho)H_R(n,m) * i * rho
    !   When ind is present only the ind component of the velocity operator
    !   is calculated otherwise all three components are calculated
    !
    !Args: (needs output from read_Hr)
    !   k: k point 
    !   ind (optional): when present, only this component of the v operator is 
    !       returned
    !
    !Returns:
    !   vk: the velocity operator at point k
    !--------------------------------------------------------------------------

    class(tb_model) :: self
    real(dp), intent(in) :: k(3)
    complex(dp), allocatable, intent(out) :: vk(:,:,:)
    integer, intent(in), optional :: ind

    integer :: n_wann
    integer :: i,j,n,m
    real(dp) :: latt_vecs(3,3)
    real(dp) :: r(3)

    n_wann = self%n_orb

    allocate(vk(n_wann,n_wann,3))
    vk(:,:,:) = 0_dp

    call read_structure(latt_vecs)

    do j=1,3
        if (present(ind)) then !if ind is present and j/=ind skip to the next cycle
            if (ind /= j) cycle
        endif
        do i = 1,self%n_Hrs
            r = self%rhos(1,i) * latt_vecs(:,1) + self%rhos(2,i) * latt_vecs(:,2) + self%rhos(3,i) * latt_vecs(:,3)
            n = self%orbs(1,i)
            m = self%orbs(2,i)
            vk(n,m,j) = vk(n,m,j) + exp(ii*dot_product(r,k)) * ii * r(j) * self%Hr(i)
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
    !   n_wann: magnitude of the Hamiltonian
    !   atom(optional): number of atom on which the spin-operator is projected
    !       if atom is not present, no projection is done
    !   projs(optional): must be present if atom is present. defines how the 
    !       wannier orbitals are localized, which defines how the projections are
    !       done. This is outputed by read_projs
    !
    !Returns:
    !   smat: all three components of the Hamiltonian
    !--------------------------------------------------------------------------
    class(tb_model) :: self
    complex(dp), allocatable, intent(out) :: smat(:,:,:)
    integer, optional, intent(in) :: atom
    type(projs_type), optional, intent(in) :: projs

    integer :: i,orb,n_wann

    n_wann = self%n_orb

    allocate(smat(n_wann,n_wann,3))
    smat(:,:,:) = 0._dp
    
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


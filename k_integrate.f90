module k_integrate
use constants
use mpi
implicit none
private
public cross_product,integrate_sum_1uc, find_recvecs,integrate_sum_1uc_para
    
contains

function cross_product(a,b) result(c)
    real(dp), intent(in) :: a(3)
    real(dp), intent(in) :: b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)

end function

subroutine find_recvecs(latt_vecs,rec_vecs,rec_vol)

    real(dp), intent(in) :: latt_vecs(3,3)
    real(dp), intent(out) :: rec_vecs(3,3)
    real(dp), intent(out) :: rec_vol

    real(dp) :: vol
    
    vol = dot_product(latt_vecs(:,1),cross_product(latt_vecs(:,2),latt_vecs(:,3)))

    rec_vecs(:,1) = 2*pi * cross_product(latt_vecs(:,2),latt_vecs(:,3)) / vol
    rec_vecs(:,2) = 2*pi * cross_product(latt_vecs(:,3),latt_vecs(:,1)) / vol
    rec_vecs(:,3) = 2*pi * cross_product(latt_vecs(:,1),latt_vecs(:,2)) / vol

    rec_vol = (2*pi)**3 / vol

end subroutine

function k2kindex(k,nk) 
    integer, intent(in) :: k(3)
    integer, intent(in) :: nk(3)
    integer :: k2kindex

    k2kindex = k(1) * nk(2)*nk(3) + k(2) * nk(3) + k(3)

end function

function kindex2k(i,nk) result(k)
    integer, intent(in) :: i
    integer, intent(in) :: nk(3)
    integer :: k(3)

    k(1) = i/(nk(2)*nk(3))
    k(2) = (i-k(1)*nk(2)*nk(3))/nk(3)
    k(3) = i-k(1)*nk(2)*nk(3)-k(2)*nk(3)

end function

subroutine distribute_k(nk,ntasks,rank,ls,lf)
    integer, intent(in) :: nk,ntasks,rank
    integer, intent(out) :: ls,lf
    
    integer :: nkp,d

    nkp = nk/ntasks
    d = nk-nkp*ntasks

    if (rank < d) then
        ls = 1+rank*nkp + rank
        lf = (1+rank)*nkp + rank + 1
    else
        ls = 1+rank*nkp + d
        lf = (1+rank)*nkp + d
    endif

end subroutine


function integrate_sum_1uc(func,nout,nk,latt_vecs)
    !--------------------------------------------------------------------------
    !Integrates an expression over k over a first primitive cell in the k space
    !
    !Args:
    !   func: the function to be integrated. func(real(dp)(3),nout)=real(dp)(nout)
    !   nout(int): number of functions to be integrated
    !   nk(integer(3)): number of k points along each direction
    !   latt_vecs(real(dp)(3,3)): Bravais lattice vectors stored as columns of
    !       latt_vecs
    !   
    !Returns:
    !   func integrated over the first primive cell
    !
    !Notes:
    !   This integral is much easier to define than the integral over the first
    !        Brillouin zone.
    !   Due to technical reasons, the function func needs as a second argument
    !       the number of integrated functions. This argument does not have to be
    !       used in the function at all, except for the definition of the function
    !       output. That is, func(k,n) have to be 1-dim array of size n.
    !   Does not contain the 1/(2*pi)^3 factor!!!
    !
    !How to use:
    !   Typically, the function func requires some more arguments apart from k.
    !   Define a wrapper like this :
    !   
    !   subroutine integrate(nk,latt_vecs,arg1,arg2,...)
    !       write(*,*) integrate_sum_1uc(func,nk,latt_vecs)
    !       
    !   contains
    !       function func(k,nout)
    !           real(dp), intent(in) :: k(3)
    !           integer, intent(in) :: nout
    !           real(dp) :: func(nout)
    !           func = func_args(k,arg1,arg2,...)
    !       end function    
    !   end subroutine
    !--------------------------------------------------------------------------


    interface
        function func(k,n)
            import
            real(dp), intent(in):: k(3)
            integer, intent(in) :: n
            real(dp) :: func(n)
        end function func
    end interface
    integer, intent(in) :: nout
    integer, intent(in) :: nk(3)
    real(dp), intent(in) :: latt_vecs(3,3)
    real(dp) :: integrate_sum_1uc(nout)

    integer :: k1,k2,k3
    real(dp) :: kf(3),k(3)
    real(dp) :: rec_vecs(3,3)
    real(dp) :: vol
    real(dp) :: vol_dk
    real(dp) :: I(nout)
    integer :: j
    real(dp) :: res(nout)

    I(:) = 0._dp

    call find_recvecs(latt_vecs,rec_vecs,vol)

    vol_dk = vol/(nk(1)*nk(2)*nk(3))

    do k1=0,nk(1)-1
        do k2=0,nk(2)-1
            do k3=0,nk(3)-1

                kf(1) = dble(k1)/dble(nk(1))
                kf(2) = dble(k2)/dble(nk(2))
                kf(3) = dble(k3)/dble(nk(3))

                k(:) = kf(1) * rec_vecs(:,1) + kf(2) * rec_vecs(:,2) + &
                    kf(3) * rec_vecs(:,3)
                
                res = func(k,nout)
                do j=1,nout
                    I(j) = I(j) + res(j)*vol_dk
                end do

            end do
        end do
    end do

    integrate_sum_1uc = I

end function        

function integrate_sum_1uc_para(func,nout,nk,latt_vecs,rank)
    !--------------------------------------------------------------------------
    !Integrates an expression over k over a first primitive cell in the k space
    !Parallel version
    !
    !Args:
    !   func: the function to be integrated. func(real(dp)(3),nout)=real(dp)(nout)
    !   nout(int): number of functions to be integrated
    !   nk(integer(3)): number of k points along each direction
    !   latt_vecs(real(dp)(3,3)): Bravais lattice vectors stored as columns of
    !       latt_vecs
    !   
    !Returns:
    !   func integrated over the first primive cell
    !
    !Notes:
    !   This integral is much easier to define than the integral over the first
    !        Brillouin zone.
    !   Due to technical reasons, the function func needs as a second argument
    !       the number of integrated functions. This argument does not have to be
    !       used in the function at all, except for the definition of the function
    !       output. That is, func(k,n) have to be 1-dim array of size n.
    !   Does not contain the 1/(2*pi)^3 factor!!!
    !
    !How to use:
    !   Typically, the function func requires some more arguments apart from k.
    !   Define a wrapper like this :
    !   
    !   subroutine integrate(nk,latt_vecs,arg1,arg2,...)
    !       write(*,*) integrate_sum_1uc(func,nk,latt_vecs)
    !       
    !   contains
    !       function func(k,nout)
    !           real(dp), intent(in) :: k(3)
    !           integer, intent(in) :: nout
    !           real(dp) :: func(nout)
    !           func = func_args(k,arg1,arg2,...)
    !       end function    
    !   end subroutine
    !--------------------------------------------------------------------------


    interface
        function func(k,n)
            import
            real(dp), intent(in):: k(3)
            integer, intent(in) :: n
            real(dp) :: func(n)
        end function func
    end interface
    integer, intent(in) :: nout
    integer, intent(in) :: nk(3)
    real(dp), intent(in) :: latt_vecs(3,3)
    integer, intent(in) :: rank
    real(dp) :: integrate_sum_1uc_para(nout)

    real(dp) :: kf(3),k(3)
    real(dp) :: rec_vecs(3,3)
    real(dp) :: vol
    real(dp) :: vol_dk
    real(dp) :: I(nout)
    integer :: j,l
    real(dp) :: res(nout)
    integer :: totnk !total number of kpoints
    integer :: ntasks
    integer :: ki(3)
    integer :: ls,lf

    !MPI vars
    integer :: ierr

    I(:) = 0._dp

    call find_recvecs(latt_vecs,rec_vecs,vol)

    vol_dk = vol/(nk(1)*nk(2)*nk(3))

    totnk = nk(1)*nk(2)*nk(3)

    call mpi_comm_size(mpi_comm_world,ntasks,ierr)

    call distribute_k(totnk,ntasks,rank,ls,lf)

    do l = ls,lf
        ki = kindex2k(l,nk)

        kf(1) = dble(ki(1))/dble(nk(1))
        kf(2) = dble(ki(2))/dble(nk(2))
        kf(3) = dble(ki(3))/dble(nk(3))

        k(:) = kf(1) * rec_vecs(:,1) + kf(2) * rec_vecs(:,2) + &
            kf(3) * rec_vecs(:,3)
        
        res = func(k,nout)
        do j=1,nout
            I(j) = I(j) + res(j)*vol_dk
        end do

    end do

    integrate_sum_1uc_para = I

end function        

end module

program mag_moms

use mpi
use constants
use input
use read, only: read_Hr_para, read_structure
use eval_k, only: create_Hk2,eigen,create_smat,mat_ele
use k_integrate, only: cross_product,integrate_sum_1uc_para

type(input_type) :: inp
type(projs_type) :: projs
integer :: n_wann
integer :: n_r
integer, allocatable :: r_deg(:)
integer, allocatable :: rpts(:,:)
complex(dp), allocatable :: Hr(:,:,:)
complex(dp), allocatable :: Hk(:,:)
real(dp) :: latt_vecs(3,3)
real(dp), allocatable :: w(:)
integer :: ierr,rank,ntasks
real(dp), allocatable :: mag_mom_p(:),mag_mom(:)
integer :: a,u

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input(inp)
call read_projs(projs)
call read_structure(latt_vecs)

call read_Hr_para(n_wann,n_r,r_deg,rpts,Hr)

allocate(mag_mom(3*inp%n_projs),mag_mom_p(3*inp%n_projs))

mag_mom_p = integrate_sum_1uc_para(func,3*inp%n_projs,inp%nk,latt_vecs,rank)
vol = dot_product(latt_vecs(:,1),cross_product(latt_vecs(:,2),latt_vecs(:,3)))

mag_mom_p = mag_mom_p * vol / (2*pi)**3 
call mpi_reduce(mag_mom_p,mag_mom,3*inp%n_projs,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

if ( rank == 0 ) then
    open(newunit=u, file='mag_mom_out', recl=300)
    do a=1,inp%n_projs
        write(u,*) 'atom', inp%atoms(a)
        write(u,*) mag_mom((a-1)*3+1:(a-1)*3+3)
        write(u,*) 
    end do
endif

call mpi_finalize(ierr)

contains 

    subroutine mag_mom_k(k,mom_k)
        real(dp), intent(in) :: k(3)
        real(dp), intent(out), allocatable :: mom_k(:)

        complex(dp),allocatable :: ele(:,:)
        complex(dp),allocatable :: smat(:,:,:)
        integer :: i,n,a,atom

        allocate(mom_k(3*inp%n_projs))
        allocate(ele(n_wann,n_wann))
        call create_Hk2(r_deg,rpts,Hr,k,Hk)
        call eigen(Hk,w)

        mom_k = 0_dp
        do i=1,3

            do a=1,inp%n_projs

                atom = inp%atoms(a)
                call create_smat(n_wann,smat,atom,projs)
                call mat_ele(smat(:,:,i),Hk,ele)

                do n=1,n_wann
                    if (w(n) < inp%Ef ) then
                        mom_k((a-1)*3+i) = mom_k((a-1)*3+i) + ele(n,n)
                    end if
                end do
            end do
        end do

    end subroutine

    function func(k,n)
        real(dp), intent(in) :: k(3)
        integer, intent(in) :: n
        real(dp) :: func(n)

        real(dp), allocatable :: mm_out(:)

        call mag_mom_k(k,mm_out)
        func = mm_out

    end function


end program

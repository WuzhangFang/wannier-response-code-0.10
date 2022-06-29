program linres
!--------------------------------------------------------------------------
!This calculates various linear response quantities using either a wannier
!tb Hamiltonian or an empirical tight-binding Hamiltonian
!--------------------------------------------------------------------------

use mpi
use constants
use input
use eval_k, only: linres_k
use k_integrate, only: cross_product,integrate_sum_1uc_para
use generic_model
use wann_model_mod
use tb_model_mod
implicit none

integer :: i,j,l,m,n,o
real(dp), allocatable :: X(:,:,:,:)
real(dp), allocatable :: Xtot(:,:,:,:)
real(dp) :: t1,t2,t3,t4
real(dp) :: times(10,2)
integer :: ierr,rank,ntasks
integer :: nktot, n_fors
integer :: n_gam, n_projs
real(dp) :: latt_vecs(3,3)
real(dp) :: vol
integer :: u
type(input_type) :: inp
type(projs_type) :: projs
integer, allocatable :: atoms(:)
class(gen_model), allocatable :: model
character(len=100), parameter :: version = '0.1.0'

character(len=100), parameter :: output_file = 'linres_out'
times = 0_dp

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
if (rank == 0) then 
    t1 = mpi_wtime()
endif
call mpi_comm_size(mpi_comm_world,ntasks,ierr)


call read_structure(latt_vecs)

call read_input(inp)

!we initialize the model chosen in the input
if (inp%model == 'wann') then
    allocate(wann_model::model)
else if (inp%model == 'tb') then
    allocate(tb_model::model)
else
    stop 'wrong model type in input'
endif
call model%init()

if (rank == 0) then 
    t2 = mpi_wtime()
endif

!if there are projections read them
if (inp%n_projs > 0) then
    call read_projs(projs)
endif

nktot = inp%nk(1)*inp%nk(2)*inp%nk(3)
!parallelization is the most efficient for commensurate number of kpoints and tasks, but 
!this doesnt matter in most cases
if (rank == 0) then
    if ((nktot/ntasks)*ntasks /= nktot) then
        write(*,*) 'Use commensurate number of kpoints and tasks for best result'
    endif
endif

n_gam = size(inp%gam)
n_fors = inp%n_even + inp%n_odd
n_projs = max(1,inp%n_projs)
allocate(X(inp%n_op1,3,n_gam,n_fors))
allocate(Xtot(inp%n_op1,3,n_gam,n_fors))

allocate(atoms(n_projs))
if (inp%n_projs == 0) then
    atoms(1) = 0
else
    atoms = inp%atoms
endif


if (rank == 0) then 
    t3 = mpi_wtime()
endif

!integrate the linear response formular over k
X = reshape(integrate_sum_1uc_para(func,inp%n_op1*3*n_gam*n_fors,inp%nk,latt_vecs,rank),[inp%n_op1,3,n_gam,n_fors])
vol = dot_product(latt_vecs(:,1),cross_product(latt_vecs(:,2),latt_vecs(:,3)))
!convert to correct units
if (inp%linres == 'cisp') then
    !Convert to units of hbar * V / m - it's a spin per unit cell per unit electric field
    X = X * vol / (2*pi)**3 * 1e-10_dp / 2  / eV
elseif (inp%linres == 'cond') then
    !convert to units 1/Omega/cm
    X = X  * (-eele) / (2*pi)**3 * 1e-8_dp / eV / 1e-16_dp
elseif (inp%linres == 'she') then
    !to get proper units we have to divide by 2 to get correct spin and convert A to cm
    !We also have to divide by 2pi^3 because thats missing in the k integral
    !then we have result in hbar/e/Omega/cm
    X = X * 10**8 /2 / (2*pi)**3
endif

call mpi_reduce(X,Xtot,inp%n_op1*3*n_gam*n_fors,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

if (rank == 0) then 
    t4 = mpi_wtime()
endif

!prints the output
if (rank == 0) then
    open(newunit=u, file=output_file, recl=300)
    write(u,*) 'Version: ', version
    write(u,*) 'total number of k-points', inp%nk(1)*inp%nk(2)*inp%nk(3)
    write(u,*) 'number of processes', ntasks
    write(u,*) 'reading input took', t2-t1, 's'
    write(u,*) 'calculation took', t4-t3, 's'
    write(u,*) 'time per 1 process and 1000 k points', (t4-t3) * ntasks / (inp%nk(1)*inp%nk(2)*inp%nk(3))
    write(u,*) 
    write(u,*) 'Detailed timings:'
    write(u,*) 'constructing the velocity operator took:', times(1,2)-times(1,1), 's'
    write(u,*) 'constructing Hk took:', times(7,2)-times(7,1), 's'
    write(u,*) 'diagonalizing the Hamiltonian took:', times(8,2)-times(8,1), 's'
    write(u,*) 'constructing the spin-operator (for all projections) took:', times(2,2)-times(2,1), 's'
    write(u,*) 'matrix elements of velocity took:', times(3,2)-times(3,1), 's'
    write(u,*) 'matrix elements of spin(for all projections) took:', times(4,2)-times(4,1), 's'
    write(u,*) 'evaluating all the linear response formulas for all gammas took:', times(5,2)-times(5,1), 's'
    write(u,*) 'linres_k total:', times(6,2)-times(6,1), 's'
    write(u,*) 
    if (inp%linres == 'cisp') then
        write(u,*) 'CISP in units of hbar*V/m'
    elseif (inp%linres == 'cond') then
        write(u,*) 'Conductivity in units of 1/(Ohm*cm)'
    elseif (inp%linres == 'she') then
        write(u,*) 'SHE in units of hbar/e/Omega/cm'
    endif

    
    do l=1,n_gam
        write(u,*) '*****************************'
        write(u,*) 'Gamma: ', inp%gam(l)
        do n=1,n_projs
            if (inp%n_projs > 0 ) then
                if (inp%atoms(n)==0) then
                    write(u,*) 'Total:'
                else
                    write(u,*) 'Projection on atom', inp%atoms(n)
                endif
            endif
            do m=1,inp%n_even
                write(u,*) 'Even part, formula: ', inp%for_even(m)
                if (inp%linres == 'cond' .or. inp%linres == 'cisp') then
                    do i =1,3
                        write(u,*) (Xtot((n-1)*3+i,j,l,m),j=1,3)
                    end do
                else if (inp%linres == 'she') then
                    do i=1,3
                        write(u,*) 'spin_component', i
                        do j=1,3
                            write(u,*) (Xtot((i-1)*3+j,o,l,m),o=1,3)
                        end do
                        write(u,*)
                    end do
                end if
                write(u,*)
            end do
            do m=1,inp%n_odd
                write(u,*) 'Odd part, formula: ', inp%for_odd(m)
                if (inp%linres == 'cond' .or. inp%linres == 'cisp') then
                    do i =1,3
                        write(u,*) (Xtot((n-1)*3+i,j,l,m+inp%n_even),j=1,3)
                    end do
                else if (inp%linres == 'she') then
                    do i=1,3
                        write(u,*) 'spin_component', i
                        do j=1,3
                            write(u,*) (Xtot((i-1)*3+j,o,l,m+inp%n_even),o=1,3)
                        end do
                        write(u,*)
                    end do
                end if
                write(u,*)
            end do
        end do
    end do
endif

!This just dumps the Xtot array without any formatting (though it is a formatted file in the fortran
!terminology)
if (rank == 0 ) then
    write(u,*) 
    write(u,*) '*****************************'
    write(u,*) 'Unformatted output:'
    do i=1,inp%n_op1
        do j=1,3
            do l=1,n_gam
                do m=1,n_fors
                    write(u,*) i,j,inp%gam(l),m,Xtot(i,j,l,m)
                end do
            end do
        end do
    end do
end if

call mpi_finalize(ierr)

contains 

    !this is a wrapper for the linres_k function so that this function can be sent to k_integrate
    function func(k,n)
        real(dp), intent(in) :: k(3)
        integer, intent(in) :: n
        real(dp) :: func(n)

        real(dp), allocatable :: Xo(:,:,:,:)

        if (inp%n_projs > 0) then 
            call linres_k(k,model,inp,Xo,projs,times)
        else
            call linres_k(k,model,inp,Xo,times=times)
        endif
        func = reshape(Xo,[inp%n_op1*3*n_gam*n_fors])

    end function


end program

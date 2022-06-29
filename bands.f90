program bands

use constants
use input
use eval_k, only: eigen
use k_integrate, only: find_recvecs
use generic_model
use wann_model_mod
use tb_model_mod

type(input_type) :: inp
type(bands_input_type) :: bands_inp
integer :: n_wann
complex(dp), allocatable :: Hk(:,:)
real(dp) :: latt_vecs(3,3)
real(dp), allocatable :: w(:)
real(dp) :: rec_vecs(3,3)
real(dp) :: rec_vol
class(gen_model), allocatable :: model

integer :: nk,i,j,n,u
real(dp), dimension (3) :: k,kf
real(dp), allocatable :: klist(:,:)

call read_input(inp)
call read_bands_input(bands_inp)
call read_structure(latt_vecs)
call find_recvecs(latt_vecs,rec_vecs,rec_vol)

if (inp%model == 'wann') then
    allocate(wann_model::model)
else if (inp%model == 'tb') then
    allocate(tb_model::model)
end if

call model%init()
n_wann = model%n_orb

nk = bands_inp%n_kpts * (bands_inp%n_sym_kpts - 1)

allocate(klist(3,nk))

do i=1,bands_inp%n_sym_kpts-1
    do j =1,bands_inp%n_kpts
        klist(:, (i-1)*bands_inp%n_kpts + j) = bands_inp%sym_kpts(:,i) + &
            (j-1)/dble(bands_inp%n_kpts) * (bands_inp%sym_kpts(:,i+1) - bands_inp%sym_kpts(:,i)) 
    end do
end do

open(newunit=u, file='bands.dat', recl=300)
do i=1,nk
    kf = klist(:,i)
    k = kf(1) * rec_vecs(:,1) + kf(2) * rec_vecs(:,2) + kf(3) * rec_vecs(:,3)
    call model%create_Hk(k,Hk)
    if ( bands_inp%print_Ham ) then
        write(u,*) ''
        do n=1,n_wann
            write(u,'(12(F8.5))') Hk(n,:)
        end do
    end if

    call eigen(Hk,w)
    do n=1,n_wann
        write(u,*) i, kf, w(n)
    end do
    write(u,*) ''
    write(u,*) ''
end do
        
end program




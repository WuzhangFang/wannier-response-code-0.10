module eval_k
use lapack95
use blas95
use mpi
use constants
use input
use generic_model
implicit none
private
public :: eigen, mat_ele,linres_k

contains

subroutine eigen(A,w)
    !--------------------------------------------------------------------------
    !Finds eigenvalues and eigenvectors of a Hermitian matrix A.
    !
    !Note:!!! Tha matrix A is overwritten!!!
    !
    !Args:
    !   A(complex(:,:)): the matrix - assumes the lower diagonal part
    !
    !Returns:
    !   A: eigenvectors - stored as the columns of A
    !   w: eigenvalues
    !--------------------------------------------------------------------------

    complex(dp),intent(inout) :: A(:,:)
    real(dp), intent(out), allocatable :: w(:)

    integer :: n
    complex(dp), allocatable :: z(:,:)

    n = size(A,1)
    allocate(w(n))
    allocate(z(n,n))

    !call heevd(A,w,'V','L')
    call heevr(A,w,z=z,uplo='L')
    A=z

end subroutine

subroutine mat_ele(mat,evecs,ele)
    !--------------------------------------------------------------------------
    !Calculates matrix elements of operator
    !
    !Args:
    !   mat: matrix of the operator
    !Returns:
    !   evecs: Matrix, which contains the eigenvectors in a row form, i.e., the
    !       rows of evecs are the eigenvectors.
    !
    !Note: Requires BLAS.
    !--------------------------------------------------------------------------

    complex(dp), intent(in) :: mat(:,:)
    complex(dp), intent(in) :: evecs(:,:)
    complex(dp), intent(out) :: ele(:,:)

    integer :: n_wann
    complex(dp), allocatable ::mat2(:,:)
    
    n_wann = size(mat,1)
    allocate(mat2(n_wann,n_wann))            
    call zgemm('N','N',n_wann,n_wann,n_wann,(1.D0,0.D0),mat(:,:),n_wann,evecs(:,:),n_wann,(0.D0,0.D0),mat2(:,:),n_wann)
    !         'C': Hermitian                   alpha       A                B                beta       C
    call zgemm('C','N',n_wann,n_wann,n_wann,(1.D0,0.D0),evecs(:,:),n_wann,mat2(:,:),n_wann,(0.D0,0.D0),ele(:,:),n_wann)

end subroutine

subroutine linres_k(k,model,inp,X,projs,times)
    !--------------------------------------------------------------------------
    !Evaluates the interband CISP at point k. When integrated over k, this
    !gives the interband cisp.
    !
    !Args;
    !   k: the k points
    !   model (class(gen_model)): can be any class that extends the gen_model
    !       class. We only need procedures for creating Hk, vk and smat from
    !       this object.
    !   inp (type(input_type)): the input type
    !       contains all the necessary information about parameters of the
    !       calculation, in particular the gamma values, Fermi level, projections
    !       and the actual quantity we are calculating
    ! 
    !Returns:
    !   X: containts the cisp tensor
    !      format: X(i,j,,Ef), where j is always the last component of the
    !       tensor (that is the component of the second operator in the
    !       linear response formula) and i labels the first operator components
    !       in case of cond: i=1,3
    !       in case of cisp: i=(n-1)*3+l, where n is the number of spin-projection
    !           and l is the first component of the cisp tensor
    !       in case of she: i=(l-1)*3 + k, where l,k are the first two components
    !           of the tensor
    !--------------------------------------------------------------------------


    real(dp), intent(in) :: k(3)
    class(gen_model) :: model
    type(input_type), intent(in) :: inp
    real(dp), intent(out), allocatable :: X(:,:,:,:)
    type(projs_type), intent(in), optional :: projs
    real(dp), optional, intent(inout) :: times(10,2)

    complex(dp), allocatable :: Hk(:,:)
    complex(dp), allocatable :: vk(:,:,:)
    complex(dp), allocatable :: smat(:,:,:)
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: vk_ele(:,:,:)
    complex(dp), allocatable :: op1_ele(:,:,:)
    complex(dp), allocatable :: svmat(:,:)
    integer :: n_wann
    integer :: n_gam
    integer :: n_projs
    integer :: n_op1
    integer, allocatable :: atoms(:)
    integer :: i,j,l,m,n
    real(dp) :: t

    !timing6<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(6,1) = times(6,1) + t
    endif
    !timing6>>>>>>>>>>>>>>>>>>>>

    n_gam = size(inp%gam)
    n_projs = max(1,inp%n_projs)
    n_wann = model%n_orb
    allocate(svmat(n_wann,n_wann))

    !timing7<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(7,1) = times(7,1) + t
    endif
    !timing7>>>>>>>>>>>>>>>>>>>>
    
    !creates the k-dependent Hamiltonian at point k
    call model%create_Hk(k,Hk)

    !timing7<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(7,2) = times(7,2) + t
    endif
    !timing7>>>>>>>>>>>>>>>>>>>>

    !timing8<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(8,1) = times(8,1) + t
    endif
    !timing8>>>>>>>>>>>>>>>>>>>>

    !finds eigenvalues and eigenfunctions of Hk
    call eigen(Hk,w)

    !timing8<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(8,2) = times(8,2) + t
    endif
    !timing8>>>>>>>>>>>>>>>>>>>>

    
    !timing1<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(1,1) = times(1,1) + t
    endif
    !timing1>>>>>>>>>>>>>>>>>>>>

    !constructs the velocity operator
    call model%create_vk(k,vk)

    !timing1<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(1,2) = times(1,2) + t
    endif
    !timing1>>>>>>>>>>>>>>>>>>>>

    allocate(vk_ele(n_wann,n_wann,3))

    !n_op1 is the number of first operators in the linear response formula
    ! conductivity: n_op1=3, cisp: n_op1=n*3 (n is number of projections), she: n_op1=9
    n_op1 = inp%n_op1

    allocate(op1_ele(n_wann,n_wann,n_op1))
    allocate(X(n_op1,3,n_gam,inp%n_even+inp%n_odd))
    X = 0_dp

    !timing2<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(2,1) = times(2,1) + t
    endif
    !timing2>>>>>>>>>>>>>>>>>>>>

    !matrix elements of the velocity operator, these are always needed
    do i=1,3
        call mat_ele(vk(:,:,i),Hk,vk_ele(:,:,i))
    end do

    !timing2<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(2,2) = times(2,2) + t
    endif
    !timing2>>>>>>>>>>>>>>>>>>>>

    allocate(atoms(n_projs))
    if (inp%n_projs == 0) then
        !if there is no projection, we set it to 0, which means total
        atoms(1) = 0
    else
        atoms = inp%atoms ! store the index of atoms
    endif

    !for cisp we construct the spin-operator for each projection
    if (inp%linres == 'cisp') then

        !a loop over the atom projections
        do n=1,n_projs

            !timing3<<<<<<<<<<<<<<<<<<<<
            if (present(times)) then
                t = mpi_wtime() 
                times(3,1) = times(3,1) + t
            endif
            !timing3>>>>>>>>>>>>>>>>>>>>

            call model%create_smat(smat,atoms(n),projs)

            !timing3<<<<<<<<<<<<<<<<<<<<
            if (present(times)) then
                t = mpi_wtime() 
                times(3,2) = times(3,2) + t
            endif
            !timing3>>>>>>>>>>>>>>>>>>>>

            !timing4<<<<<<<<<<<<<<<<<<<<
            if (present(times)) then
                t = mpi_wtime() 
                times(4,1) = times(4,1) + t
            endif
            !timing4>>>>>>>>>>>>>>>>>>>>

            !loops over the components of the cisp tensor
            do i =1,3
                call mat_ele(smat(:,:,i),Hk,op1_ele(:,:,(n-1)*3+i))
            end do

            !timing4<<<<<<<<<<<<<<<<<<<<
            if (present(times)) then
                t = mpi_wtime() 
                times(4,2) = times(4,2) + t
            endif
            !timing4>>>>>>>>>>>>>>>>>>>>

        end do

    else if (inp%linres == 'cond') then
        !for conductivit the first operator is just the velocity operator
        op1_ele = vk_ele
    else if (inp%linres == 'she') then
        !for SHE the first operator is s*v

        !timing3<<<<<<<<<<<<<<<<<<<<
        if (present(times)) then
            t = mpi_wtime() 
            times(3,1) = times(3,1) + t
        endif
        !timing3>>>>>>>>>>>>>>>>>>>>

        call model%create_smat(smat)

        !timing3<<<<<<<<<<<<<<<<<<<<
        if (present(times)) then
            t = mpi_wtime() 
            times(3,2) = times(3,2) + t
        endif
        !timing3>>>>>>>>>>>>>>>>>>>>

        do i=1,3
            do j=1,3
                !calculates 1/2 of anticommutator of s and v, 0.5*(s*v+v*s)
                ! C := alpha*op(A)*op(B) + beta*C
                !           A         B          C             alpha        beta
                call gemm(vk(:,:,j),smat(:,:,i),svmat,alpha=dcmplx(0.5_dp),beta=dcmplx(0))
                call gemm(smat(:,:,i),vk(:,:,j),svmat,alpha=dcmplx(0.5_dp),beta=dcmplx(1))
                call mat_ele(svmat,Hk,op1_ele(:,:,3*(i-1)+j))
            end do
        end do
    else
        write(*,*) inp%linres, 'is a wrong switch'
        stop 'wrong switch'
    endif

    !timing5<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(5,1) = times(5,1) + t
    endif
    !timing5>>>>>>>>>>>>>>>>>>>>

    !a loop over all the even formulas
    do m=1,inp%n_even
        !for this formula, the gamma dependence is 1/Gamma so we don't have to calculate for
        !all gammas
        if (inp%for_even(m) == 1) then
            X(:,:,1,m) = even_eval(inp%gam(1),inp%Ef,w,op1_ele,vk_ele,inp%for_even(m),inp%eps)
            do l=2,n_gam
                X(:,:,l,m) = X(:,:,1,m) * inp%gam(1) / inp%gam(l)
            end do

        !for other formulas we loop over gamma
        else
            do l=1,n_gam
                X(:,:,l,m) = even_eval(inp%gam(l),inp%Ef,w,op1_ele,vk_ele,inp%for_even(m),inp%eps)
            end do
        endif
    end do

    !a loop over all the odd formulas
    do m=1,inp%n_odd
        !a loop over gamma values
        do l=1,n_gam
            X(:,:,l,m+inp%n_even) = odd_eval(inp%gam(l),inp%Ef,w,op1_ele,vk_ele, inp%for_odd(m))
        end do
    end do

    !timing5<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(5,2) = times(5,2) + t
    endif
    !timing5>>>>>>>>>>>>>>>>>>>>

    !timing6<<<<<<<<<<<<<<<<<<<<
    if (present(times)) then
        t = mpi_wtime() 
        times(6,2) = times(6,2) + t
    endif
    !timing6>>>>>>>>>>>>>>>>>>>>

end subroutine

function even_eval(gam,Ef,w,mat1,mat2,formula,eps) result(x)
    !--------------------------------------------------------------------------
    !Evaluates the interband linear response formula at point k.
    !--------------------------------------------------------------------------

    real(dp), intent(in) :: gam
    real(dp), intent(in) :: Ef
    real(dp), intent(in) :: w(:)
    complex(dp), intent(in) :: mat1(:,:,:)
    complex(dp), intent(in) :: mat2(:,:,:)
    integer, intent(in) :: formula
    real(dp), intent(in), optional :: eps

    real(dp), allocatable :: x(:,:)

    integer :: i,j,n,m,n_wann
    integer :: nmat1,nmat2
    real(dp) :: f
    complex(dp), allocatable :: vec1(:),vec2(:),vec3(:),vec4(:)
    complex(dp), allocatable :: mat3(:,:)
    real(dp) :: t(5,2),tt
    integer, allocatable :: order(:)

    t(1,1) = mpi_wtime()
    n_wann = size(w)
    allocate(vec1(n_wann**2),vec2(n_wann**2),vec3(n_wann**2),vec4(n_wann**2))
    allocate(mat3(n_wann,n_wann))

    nmat1 = size(mat1,3)
    nmat2 = size(mat2,3)

    if (.not. allocated(x)) then
        allocate(x(nmat1,nmat2))
    endif
    x = 0._dp

    if (formula == 1) then
        if (.not. present(eps)) stop 'you must provide epsilon'
        do n=1,n_wann
            do i=1,nmat1
                do j=1,nmat2
                    x(i,j) = x(i,j) + mat1(n,n,i) * mat2(n,n,j) * eps / ( eps**2 + (w(n)-Ef)**2 ) / pi
                end do
            end do
        end do
        x = x * (-eele) * hbar / (2*gam)
    elseif (formula == 2) then
        do n=1,n_wann
            do m=1,n_wann
                f = gam**2 / ( ((Ef-w(n))**2 + gam**2) * ((Ef-w(m))**2 + gam**2) )
                do i=1,nmat1
                    do j=1,nmat2
                        x(i,j) = x(i,j) + dble(mat1(n,m,i)*mat2(m,n,j)) * f ! dble(a) converts a to a real number
                    end do
                end do
            end do
        end do
        x = x * (-eele)*hbar/pi

    !just some testing switch
    elseif (formula == 21) then

        tt = mpi_wtime()
        t(2,1) = t(2,1) + tt
        do n=1,n_wann
            do m=1,n_wann
                mat3(n,m) =  gam**2 / ( ((Ef-w(n))**2 + gam**2) * ((Ef-w(m))**2 + gam**2) )
            end do
        end do

        allocate(order(n_wann**2))
        do n=1,n_wann**2
            order=n_wann**2+1-n
        end do
        vec3 = reshape(mat3,[n_wann**2])
        tt = mpi_wtime()
        t(2,2) = t(2,2) + tt

        do i=1,nmat1
            do j=1,nmat2

                tt = mpi_wtime()
                t(3,1) = t(3,1) + tt
                vec1 = reshape(mat1(:,:,i),[n_wann**2])
                vec2 = reshape(mat2(:,:,j),[n_wann**2],order=order)
                tt = mpi_wtime()
                t(3,2) = t(3,2) + tt

                tt = mpi_wtime()
                t(4,1) = t(4,1) + tt
                call vzmul(n_wann**2,vec2,vec3,vec4)
                x(i,j) = dotu(vec1,vec4)
                tt = mpi_wtime()
                t(4,2) = t(4,2) + tt

            end do
        end do

        write(*,*) 'mat3', t(2,2)-t(2,1)
        write(*,*) 'reshape', t(3,2)-t(3,1)
        write(*,*) 'eval', t(4,2)-t(4,1)
        x = x * eele*hbar/pi

    !just some testing switch
    elseif (formula == 22) then

!        t(2,1) = mpi_wtime()
        do n=1,n_wann
            do m=1,n_wann
                mat3(n,m) =  gam**2 / ( ((Ef-w(n))**2 + gam**2) * ((Ef-w(m))**2 + gam**2) )
            end do
        end do

!        t(2,2) = mpi_wtime()

    !just some testing switch
        t(3,1) = mpi_wtime()
        write(*,*) nmat1,nmat2
        do i=1,nmat1
            do j=1,nmat2

                do n=1,n_wann
                    do m=1,n_wann
                        x(i,j) = x(i,j) + mat1(n,m,i)*mat2(m,n,j)*mat3(n,m)
                    end do
                end do

            end do
        end do
        t(3,2) = mpi_wtime()

        x = x * eele*hbar/pi
    elseif (formula /= 0) then
        stop 'wrong intraband switch'
    endif


end function


function odd_eval(gam,Ef,w,mat1,mat2,formula) result(x)
    !--------------------------------------------------------------------------
    !Evaluates the interband linear response formula at point k.
    !
    !Description:
    !   The expression evaluated is from Phys. Rev. B 91, 134402 (2015), Eq. 4:
    !   2e\hbar \sum_{n\neq m} Im(mat1_nm * mat2_mn) (gam**2-(E_kn-E_km)**2)
    !       / ((E_kn-E_km)**2+gam**2) ; for n, m such that E_kn < Ef and E_km> Ef
    !   This is equivalent to: 

    !   e\hbar \sum_{n\neq m} Im(mat1_nm * mat2_nm)(gam**2-(E_kn-E_km)**2)
    !       / ((E_kn-E_km)**2+gam**2) * (f_kn-f_km)
    !   and also to:
    !   2e\hbar \sum_{n\neq m} Im(mat1_nm * mat2_nm) (gam**2-(E_kn-E_km)**2)
    !       / ((E_kn-E_km)**2+gam**2) * f_kn
    !
    !Args:
    !   gam: value of Gamma
    !   Ef: the Fermi level
    !   w: the eigenvalues
    !   mat1: matrix elements of operator 1
    !   mat2: matrix elements of operator 2
    !
    !Returns:
    !   x: the result
    !--------------------------------------------------------------------------

    real(dp), intent(in) :: gam
    real(dp), intent(in) :: Ef
    real(dp), intent(in) :: w(:)
    complex(dp), intent(in) :: mat1(:,:,:)
    complex(dp), intent(in) :: mat2(:,:,:)
    integer, intent(in) :: formula

    real(dp), allocatable :: x(:,:)

    integer :: n_wann
    integer :: n,m,i,j
    integer :: nmat1,nmat2
    real(dp) :: f,f2,f3

    n_wann = size(w)

    nmat1 = size(mat1,3)
    nmat2 = size(mat2,3)

    allocate(x(nmat1,nmat2))

    x = 0._dp

    if (formula == 1) then
        do n=1,n_wann
            do m=1,n_wann
                if ( n /= m) then
                    if ((w(n) < Ef .and. w(m) > Ef)) then
                        f =  (-( w(n) - w(m) ) ** 2.D0 - gam ** 2.D0 ) / ( ( ( w(n) - w(m) ) ** 2.D0 + gam ** 2.D0 ) ** 2.D0 )
                        do i=1,nmat1
                            do j=1,nmat2
                                x(i,j) = x(i,j) + aimag(mat1(n,m,i)*mat2(m,n,j)) * f 
                            end do
                        end do
                    end if
                endif
            end do
        end do
        x = x * 2
        x = x * eele * hbar
    elseif (formula == 2) then
        do n=1,n_wann
            do m=1,n_wann
                if ( n /= m) then
                    f =  gam * (w(m)-w(n)) / ( ((Ef-w(n))**2 + gam**2) * ((Ef-w(m))**2 + gam**2) )
                    f2 = 2 * gam / ( (w(n)-w(m)) * ((Ef-w(m))**2 + gam**2) )
                    f3 = 2 / (w(n)-w(m))**2 * aimag(log((w(m)-Ef-ii*gam) / (w(n)-Ef-ii*gam)))
                    do i=1,nmat1
                        do j=1,nmat2
                            x(i,j) = x(i,j) + aimag(mat1(n,m,i)*mat2(m,n,j)) * (f + f2 + f3) 
                        end do
                    end do 
                endif
            end do
        end do
        x = x * eele * hbar / (2*pi)
    elseif (formula /= 0) then
        stop 'wrong interband switch'
    endif
        
end function


function inter_eval2(gam,Ef,w,mat1,mat2) result(x)
    !--------------------------------------------------------------------------
    !Same as inter_eval, except it evaluates the formula:
    !   2e\hbar \sum_{n\neq m} Im(mat1_nm * mat2_nm) (gam**2-(E_kn-E_km)**2)
    !       / ((E_kn-E_km)**2+gam**2) * f_kn
    !--------------------------------------------------------------------------

    real(dp), intent(in) :: gam
    real(dp), intent(in) :: Ef
    real(dp), intent(in) :: w(:)
    complex(dp), intent(in) :: mat1(:,:)
    complex(dp), intent(in) :: mat2(:,:)

    real(dp) :: x

    integer :: n_wann
    integer :: n,m
    real(dp) :: f

    n_wann = size(w)

    x = 0._dp
    do n=1,n_wann
        do m=1,n_wann
            if ( n /= m) then

                if (w(n) < Ef) then
                    f =  (( w(n) - w(m) ) ** 2.D0 - gam ** 2.D0 ) / ( ( w(n) - w(m) ) ** 2.D0 + gam ** 2.D0 ) ** 2.D0
                    x = x + aimag(mat1(n,m)*mat2(m,n)) * f 
                end if

            endif
        end do
    end do

    x = x * 2
    x = x * eele * hbar

end function



end module


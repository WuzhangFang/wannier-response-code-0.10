% This is function to reproduce the Fortran subroutine below:
% % % subroutine mat_ele(mat,evecs,ele)
% % %     !--------------------------------------------------------------------------
% % %     !Calculates matrix elements of operator
% % %     !
% % %     !Args:
% % %     !   mat: matrix of the operator
% % %     !Returns:
% % %     !   evecs: Matrix, which contains the eigenvectors in a row form, i.e., the
% % %     !       rows of evecs are the eigenvectors.
% % %     !
% % %     !Note: Requires BLAS.
% % %     !--------------------------------------------------------------------------
% % % 
% % %     complex(dp), intent(in) :: mat(:,:)
% % %     complex(dp), intent(in) :: evecs(:,:)
% % %     complex(dp), intent(out) :: ele(:,:)
% % % 
% % %     integer :: n_wann
% % %     complex(dp), allocatable ::mat2(:,:)
% % %     
% % %     n_wann = size(mat,1)
% % %     allocate(mat2(n_wann,n_wann))
% % %     
% % %     zgemm calculates C:=alpha*op(A)*op(B)+beta*C
% % %                                              alpha       A                B                  beta         C                                                                                  
% % %     call zgemm('N','N',n_wann,n_wann,n_wann,(1.D0,0.D0),mat(:,:),n_wann,evecs(:,:),n_wann,(0.D0,0.D0),mat2(:,:),n_wann)
% % %            
% % %                 Hermitian                    alpha       A                B                  beta         C    
% % %     call zgemm('C','N',n_wann,n_wann,n_wann,(1.D0,0.D0),evecs(:,:),n_wann,mat2(:,:),n_wann,(0.D0,0.D0),ele(:,:),n_wann)
% % %     <evecs|mat|evevs>
% % % end subroutine

load('tbhr.mat');
n_band = 6;

% lattice vectors
latt_vecs = [1 0 0;
             0 1 0;
             0 0 1;];
         
k = [1 0 0];

Hk = create_Hk(k, n_band, tbhr, latt_vecs);
[Hk, w] = eig(Hk, 'vector');

vk = create_vk(k, n_band, tbhr, latt_vecs);
vk_ele = zeros(size(vk));
for i=1:3
    vk_ele(:,:,i)=Hk'*vk(:,:,i)*Hk;
end %for

m_111 = vk_ele(1,1,1); % same with Fortran (604513754757010.,3.125000000000000E-002)
m_661 = vk_ele(6,6,1); % same with Fortran (-1.687602705981830E+015,0.000000000000000E+000)

m_311 = vk_ele(3,1,1); % (468307482206482.,693419917628407.)
f_311 = complex(468307482206482,693419917628407);
norm(m_311)
norm(f_311)



function ele = mat_ele(mat,evecs,ele)
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

end
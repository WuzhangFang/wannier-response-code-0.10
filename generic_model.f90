module generic_model
integer, parameter,private :: dp=kind(0.d0)

type, abstract :: gen_model
    integer :: n_orb
    contains
        procedure(init_Int), deferred :: init
        procedure(create_smat_Int), deferred :: create_smat
        procedure(create_Hk_Int), deferred :: create_Hk
        procedure(create_vk_Int), deferred :: create_vk
end type

abstract interface 

    subroutine init_Int(self)
        import gen_model
        class(gen_model) :: self
    end subroutine

    subroutine create_smat_Int(self,smat,atom,projs)
        use input, only: projs_type
        import gen_model
        import dp
        class(gen_model) :: self
        complex(dp), allocatable, intent(out) :: smat(:,:,:)
        integer, optional, intent(in) :: atom
        type(projs_type), optional, intent(in) :: projs
    end subroutine

    subroutine create_Hk_Int(self,k,Hk)
        import gen_model
        import dp
        class(gen_model) :: self
        real(dp), intent(in) :: k(3)
        complex(dp), allocatable, intent(out) :: Hk(:,:)
    end subroutine

    subroutine create_vk_Int(self,k,vk,ind)
        import gen_model
        import dp
        class(gen_model) :: self
        real(dp), intent(in) :: k(3)
        complex(dp), allocatable, intent(out) :: vk(:,:,:)
        integer, intent(in), optional :: ind
    end subroutine

end interface

end module

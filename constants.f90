module constants
implicit none

!constants, units, etc
integer, parameter :: dp=kind(0.d0)
real(dp),  parameter :: pi  = 4 * atan (1.0_dp)
complex(dp), parameter :: ii = dcmplx(0_dp,1_dp)
real(dp), parameter :: hbar = 6.582119514e-16_dp ![eVs] Reduced Planck's constants
real(dp), parameter :: eele = 1.6021766208e-19_dp ![C] Elementary charge
real(dp), parameter :: eV = 1.602176565e-19_dp !eV

!some parameters
character(len=100), parameter :: input_file = 'input'
character(len=100), parameter :: Hr_file = 'wannier90_hr.dat'
character(len=100), parameter :: struct_file = 'POSCAR'
character(len=100), parameter :: projs_file = 'wann_projs'
character(len=100), parameter :: bands_input_file = 'input_bands'


end module

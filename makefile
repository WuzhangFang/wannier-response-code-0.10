SHELL = /bin/sh

link =  ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
compile =  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
objects = constants.o input.o generic_model.o wann_model.o tb_model.o eval_k.o k_integrate.o

linres.x: $(objects) linres.o
	mpif90 -o linres.x $(objects) linres.o -warn all -check all -traceback $(link)

bands.x: $(objects) bands.o
	mpif90 -o linres.x $(objects) bands.o -warn all -check all -traceback $(link)

mag_moms.x: $(objects) mag_moms.o
	mpif90 -o linres.x $(objects) mag_moms.o -warn all -check all -traceback $(link)

%.o: %.f90
	mpif90 -c $< $(compile)

.PHONY: clean
clean:
	rm *.o *.mod
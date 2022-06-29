# Wannier linear response
[Wannier linear response](https://bitbucket.org/zeleznyj/wannier-linear-response/wiki/Home) is a code for calculating linear response properties based on a tight-binding Hamiltonian from Wannier90. Here I try to review the first version 0.10 to learn about how the linear response formulas are implemented. Further I use Matlab to rewrite the main program to be familiar with the work flow and subroutines.

## linres.f90    

1. read the input: call read_input(inp), call read_structure(latt_vecs), call read_projs(projs)
2. choose the tight binding model: allocate(wann_model::model)
3. choose the quantity to calculate: call linres_k(k,model,inp,Xo,projs,times), then call integrate_sum_1uc_para() to integrate linres_k() over all the kpoints.

## linres_k.f90

1. creates the k-dependent Hamiltonian at point k: call model%create_Hk(k,Hk)

2. finds eigenvalues and eigenfunctions of Hk: call eigen(Hk,w). after calling, Hk becomes columns of eigenvectors.

3. constructs the velocity operator: call model%create_vk(k,vk)

4. allocates the matrices: allocate(vk_ele(n_wann,n_wann,3)), allocate(op1_ele(n_wann,n_wann,n_op1)), allocate(X(n_op1,3,n_gam,inp%n_even+inp%n_odd))

5. matrix elements of the velocity operator: call mat_ele(vk(:,:,i),Hk,vk_ele(:,:,i)). calculate matrix elements vk_ele of operator vk in the eigenvectors of Hk.

6. matrix elements of cisp tensor: call model%create_smat(smat,atoms(n),projs), call mat_ele(smat(:,:,i),Hk,op1_ele(:,:,(n-1)*3+i)). when calculating cisp, there are projections of selected orbitals.

7. matrix elements of conductivity tensor: op1_ele = vk_ele

8. matrix elements of spin current tensor: call model%create_smat(smat), call mat_ele(svmat,Hk,op1_ele(:,:,3*(i-1)+j))

9. loop over the even formulas

10. loop over the odd formulas

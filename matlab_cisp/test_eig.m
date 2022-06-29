load('tbhr.mat');
n_band = 6;

% lattice vectors
latt_vecs = [1 0 0;
             0 1 0;
             0 0 1;];
         
k = [1 0 0];

A = create_Hk(k, n_band, tbhr, latt_vecs);
[V, D] = eig(A, 'vector');

% from Fortran code the first eigenvector
a = zeros(6,1);
a(1,1)=complex(-0.209093977813657,0.000000000000000E+000);
a(2,1)=complex(-9.057960320595454E-002,-0.359206115993816);
a(3,1)=complex(0.305523425262524,0.382404493501473);
a(4,1)=complex(-0.570327366554888,-0.183184258935583);
a(5,1)=complex(0.156650798436565,0.182007391120729);
a(6,1)=complex(0.403559690378219,-1.057759843837291E-002);
b = V(:,1);
norm(b(1))
norm(a(1))
% test result
% eigen values are same with Fortran code.
% eigen vectors are not same but may be not unique.
% the norm is almost same.
% The routine is heevr from Lapack.
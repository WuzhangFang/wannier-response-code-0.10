load('tbhr.mat');
n_band = 6;

Ef = 0; 
gamma = 1e-3;

% lattice vectors
latt_vecs = [1 0 0;
             0 1 0;
             0 0 1;];
   
% test 1         
k1 = [1 0 0];
X1 = linres_k(k1,n_band,tbhr,latt_vecs,gamma,Ef);
% same with Fortran
%  -6.011429814031862E-012  0.000000000000000E+000  0.000000000000000E+000
%   0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
%   0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000

% test2
k2 = [0 1 0];
X2 = linres_k(k2,n_band,tbhr,latt_vecs,gamma,Ef);
% same with Fortran
%   0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
%   0.000000000000000E+000 -6.011429814031849E-012  0.000000000000000E+000
%   0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000

% test3
k3 = [0.5 0.5 0.5];
X3 = linres_k(k3,n_band,tbhr,latt_vecs,gamma,Ef);
% same with Fortran
%  -8.403907654708261E-013 -3.732543794501437E-013 -3.732543794501432E-013
%  -3.732543794501437E-013 -8.403907654708251E-013 -3.732543794501425E-013
%  -3.732543794501432E-013 -3.732543794501425E-013 -8.403907654708245E-013

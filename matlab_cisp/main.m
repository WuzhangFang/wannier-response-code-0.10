% Here we use a wannier Hamiltonian (tb_hr.dat) to learn how to 
% implement the Kubo formula (linres_k.m) to calculate the current-induced
% spin polarization (cisp). The Kubo formula is very similar to the one 
% to calculate the conductivity. All we need to change is to insert the
% correct operator. 
% Since the reading of Hamiltonian of Wannier format takes a lot of time 
% and memory in matlab, here we use the simple tight-binding Hamiltonian.
% Since the value is too small, I tested the linres_k.m only. It agrees
% well with the fortran result.
% @test:
% reference value for even part:
% nk = [10 10 10]
%  CISP in units of hbar*V/m
%  *****************************
%  Gamma:   1.000000000000000E-003
%  Even part, formula:            2
%  -2.181200402936098E-028 -9.371667554025620E-028 -7.440949740562756E-028
%  -2.632823271175127E-028  2.141757357675047E-028  6.840912614603759E-028
%   4.837689501267855E-028  7.502067208651822E-028  6.154360936151072E-029

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hamiltonian matrix is store in the tbhr.mat.
% Each coloumn is: R1, R2, R3, i, j, h_real, h_imaginary

% parameters: here we just input the parameters manually

load('tbhr.mat');
n_band = 6;

% lattice vectors
latt_vecs = [1 0 0;
             0 1 0;
             0 0 1];
% get the reciprocal vectors         
[rec_vecs, rec_vol, vol] = find_recvecs(latt_vecs);         
         
nk = [2 1 1];
Ef = 0; 
gamma = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start the calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum of linres_k(k) over the nk kpoints in the BZ
vol_dk = rec_vol/(nk(1)*nk(2)*nk(3));
I = zeros(9,1);

% loop over the kpoints
for k1=0:nk(1)-1
    for k2=0:nk(2)-1
        for k3=0:nk(3)-1
            k = (k1 * rec_vecs(:,1)) / nk(1) + ...
                (k2 * rec_vecs(:,2)) / nk(2) + ...
                (k3 * rec_vecs(:,3)) / nk(3);
          
            res = linres_k(k,n_band,tbhr,latt_vecs,gamma,Ef);
            for j=1:9
                I(j) = I(j) + res(j)*vol_dk;                
            end
        end %for
    end %for
end %for

X = reshape(I,[3,3]);
% add the unit hbar * V / m

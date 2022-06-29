% Here we use a Hamiltonian of s-d model (tb_hr.dat) to learn how to 
% implement the Kubo formula (linres_k.m) to calculate the conductivity. 

% @test:
% reference value for even part:
% nk = [10 10 10], gamma = 1e-3, conductivity = 80.8082615970544 
% in unit of 1/(Ohm*cm)
% test successful on 10/20/2020

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
             0 0 1;];

% get the reciprocal vectors         
[rec_vecs, rec_vol, ~] = find_recvecs(latt_vecs);         
         
nk = [10 10 10];
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
% add the unit 1/(Ohm*cm)
X = X  * -403144.194455547;

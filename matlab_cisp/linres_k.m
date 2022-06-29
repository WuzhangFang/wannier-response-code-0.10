function X = linres_k(k,n_band,tbhr,latt_vecs,gamma,Ef)
% @paramters: k, 1x3 vector
% @constants: 
hbar = 6.582119514e-16; % 
eele = 1.602176e-19; % elementary charge
% for conductivity
%X = zeros(3,3);
% for cisp
X = zeros(3,3);

% create Hk, vk and spin operator
Hk=create_Hk(k, n_band, tbhr, latt_vecs);
vk=create_vk(k, n_band, tbhr, latt_vecs);
smat = create_smat(n_band);

% get the eigenvectors and eigenvalues of H(k)
[Hk,w] = eig(Hk,'vector');

% get the matrix elements of vk operator
vk_ele = zeros(size(vk));
for i=1:3
    vk_ele(:,:,i)=Hk'*vk(:,:,i)*Hk;
end %for

% get the matrix elements of spin operator
smat_ele = zeros(size(smat));
for i=1:3
    smat_ele(:,:,i)=Hk'*smat(:,:,i)*Hk;
end %for

% formula 2
for n=1:n_band
    for m=1:n_band
        f = gamma^2 / (((Ef-w(n))^2 + gamma^2) * ((Ef-w(m))^2 + gamma^2));
        for i=1:3
            for j=1:3
                X(i,j) = X(i,j) + real(smat_ele(n,m,i)*vk_ele(m,n,j)) * f;
            end
        end
    end
end
X = X * (-eele) * hbar / pi;
end
function smat = create_smat(n_band)
% get the spin operator
% Assumes that first half othe basis functions are spin-up and the other half
% are spin-down.
% dimensionless, i.e. S=sigma. Divide by 2 to get spin in the units of hbar
% sigmax=[0 1; sigmay=[0 -i; sigmaz=[1 0; 
%         1 0]         i  0]         0 -1]
% @test
% compared to Fortran code

smat = zeros(n_band, n_band, 3);

% sigma_x
for i=1:n_band/2
    smat(i,n_band/2+i,1) = 1;
    smat(n_band/2+i,i,1) = 1;
end

% sigma_y
for i=1:n_band/2
    smat(i,n_band/2+i,2) = -1j;
    smat(n_band/2+i,i,2) = 1j;
end

% sigma_z
for i=1:n_band/2
    smat(i,i,3) = 1;
    smat(n_band/2+i,n_band/2+i,3) = -1;
end

end % function
function Hk = create_Hk(k, n_band, tbhr, latt_vecs)
% get the Hamiltonian at k
% compared to Fortran code

n_Hrs = size(tbhr, 1);
Hk = zeros(n_band, n_band);

for i=1:n_Hrs
    r = tbhr(i,1:3) * latt_vecs';
    n = tbhr(i,4);
    m = tbhr(i,5);
    Hr = complex(tbhr(i,6), tbhr(i,7));
    Hk(n,m) = Hk(n,m) + exp(1j*dot(k, r)) * Hr;
end % for


end % function
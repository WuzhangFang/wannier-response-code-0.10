function vk = create_vk(k, n_band, tbhr, latt_vecs)
% get the velocity operator at k
% a derivative of Hk
% @test
% vk1 = create_vk([1 0 0], n_band, tbhr, latt_vecs)
% vk1(1,2,1)=-1.4568e15
% vk1(2,3,1)=-1.4568e15
% vk1(6,6,3)=0
% compared to Fortran code

hbar = 6.582119514e-16;
n_Hrs = size(tbhr, 1);
vk = zeros(n_band, n_band, 3);

for j=1:3
    for i=1:n_Hrs
        r = tbhr(i,1:3) * latt_vecs';
        n = tbhr(i,4);
        m = tbhr(i,5);
        Hr = complex(tbhr(i,6), tbhr(i,7));
        vk(n,m,j) = vk(n,m,j) + exp(1j*dot(k, r)) * Hr * 1j * r(j);
    end % for
end %for

vk = vk / hbar;

end % function
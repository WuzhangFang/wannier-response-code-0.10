load('tbhr.mat');
n_band = 6;

% lattice vectors
latt_vecs = [1 0 0;
             0 1 0;
             0 0 1;];
         
k = [1 0 0];

Hk=create_Hk(k, n_band, tbhr, latt_vecs);
vk=create_vk(k, n_band, tbhr, latt_vecs);


vk(1,2,1)
vk(2,3,1)
vk(6,6,3)
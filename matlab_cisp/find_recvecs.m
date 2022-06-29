function [rec_vecs, rec_vol, vol] = find_recvecs(latt_vecs)
    % find the reciprocal lattice vectors, volume in reciprocal and real
    % space
    % @parameters: latt_vecs: lattice vectors
    vol = dot(latt_vecs(:,1), cross(latt_vecs(:,2), latt_vecs(:,3)));
    
    rec_vecs = zeros(size(latt_vecs));
    rec_vecs(:,1) = 2*pi * cross(latt_vecs(:,2), latt_vecs(:,3)) / vol;
    rec_vecs(:,2) = 2*pi * cross(latt_vecs(:,3), latt_vecs(:,1)) / vol;
    rec_vecs(:,3) = 2*pi * cross(latt_vecs(:,1), latt_vecs(:,2)) / vol;
    
    
    rec_vol = (2*pi)^3 / vol;

end
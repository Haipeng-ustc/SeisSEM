function [ W ] = WMatrix(NSPEC, NGLL, rho, vs)
%WMATRIX Summary of this function goes here
%   Local contributions of stiffness matrix (mu*wgll2): the material
%   properties and the geometry will go here

    mu = zeros(NGLL,NGLL);
    [xgll,wgll,H, HW] = GetGLL(NGLL);
    wgll2 = wgll * wgll' ;
    W     = zeros(NGLL, NGLL, NSPEC);
    for ispec = 1 : NSPEC
%         mu(:,:) = rho(:,:,ispec).*vs(:,:,ispec).^2;
        % add here the properties of heterogeneous medium
        W(:,:,ispec) = wgll2.*vs(:,:,ispec).^2;
    end


end


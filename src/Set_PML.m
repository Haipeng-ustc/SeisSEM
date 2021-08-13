function [pml_dx, pml_dz, pml_dxx, pml_dzz] = Set_PML( vp, xmin, xmax, zmin, zmax, NSPEC, NGLOB, NGLLX, NGLLZ, pml_x_thick, pml_z_thick, free_surface, coord, ibool)

USE_PML_X_Left  = true;
USE_PML_X_Right = true;
USE_PML_Z_Bottom = true;

if ( free_surface == 1 )
    USE_PML_Z_Top = false;
else
    USE_PML_Z_Top = true;
end 

Rcoef = 0.0001;
xmin_pml = xmin + pml_x_thick;
xmax_pml = xmax - pml_x_thick;
zmin_pml = zmin + pml_z_thick;
zmax_pml = zmax - pml_z_thick;

vp_max = max(vp(:));
d0_x = 3.0 * vp_max * log10(1.0/Rcoef) / (2.0 * pml_x_thick);
d0_z = 3.0 * vp_max * log10(1.0/Rcoef) / (2.0 * pml_z_thick);

% pml_dx   = zeros(NGLLX, NGLLZ, NSPEC);
% pml_dz   = zeros(NGLLX, NGLLZ, NSPEC);
% pml_dxx  = zeros(NGLLX, NGLLZ, NSPEC);
% pml_dzz  = zeros(NGLLX, NGLLZ, NSPEC);

pml_dx   = zeros(NGLOB,1);
pml_dz   = zeros(NGLOB,1);
pml_dxx  = zeros(NGLOB,1);
pml_dzz  = zeros(NGLOB,1);

%%%%%%%%%% x direction %%%%%%%%
for i = 1 : NGLOB

    x = coord(1,i);
    z = coord(2,i);

    if ( (x <= xmin_pml) && USE_PML_X_Left )
        pml_dx_temp  =   d0_x * (   (xmin_pml-x) / pml_x_thick)^2;
        pml_dxx_temp = - d0_x * 2 * (xmin_pml-x) / pml_x_thick ^2; % Be careful when caculate the derivative, expecially the direction

    elseif ( (x >= xmax_pml) && USE_PML_X_Right )
        pml_dx_temp  =   d0_x * (   (x-xmax_pml) / pml_x_thick)^2;
        pml_dxx_temp =   d0_x * 2 * (x-xmax_pml) / pml_x_thick ^2;
    else 
        pml_dx_temp  = 0.0;
        pml_dxx_temp = 0.0;
    end 

    %%%%%%%%%% z direction %%%%%%%%
    if ( (z <= zmin_pml) && USE_PML_Z_Bottom )
        pml_dz_temp  =   d0_z * (   (zmin_pml-z) / pml_z_thick)^2;
        pml_dzz_temp = - d0_z * 2 * (zmin_pml-z) / pml_z_thick ^2; % Be careful when caculate the derivative, expecially the direction
    elseif ( (z >= zmax_pml) && USE_PML_Z_Top )
        pml_dz_temp  =   d0_z * (   (z-zmax_pml) / pml_z_thick)^2;
        pml_dzz_temp =   d0_z * 2 * (z-zmax_pml) / pml_z_thick ^2;
    else 
        pml_dz_temp  = 0.0;
        pml_dzz_temp = 0.0;
    end 
    
    % for the global matrix
%     index = (ibool==i);
    pml_dx(i)  = pml_dx_temp;
    pml_dz(i)  = pml_dz_temp;
    pml_dxx(i) = pml_dxx_temp;
    pml_dzz(i) = pml_dzz_temp;

end

% dx = dx';
% dz = dz';
% dxx = dxx';
% dzz = dzz';

end 
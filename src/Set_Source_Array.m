function [sourcearray, ispec_source] = Set_Source_Array(src_node, src_type, force_angle, Moment_Tensor, ibool, xix, xiz, gammax, gammaz, NDIM, NGLLX, NGLLZ)

[xigll,wxgll,hprime_xx, hprimewgll_xx] = GetGLL(NGLLX);
[zigll,wzgll,hprime_zz, hprimewgll_zz] = GetGLL(NGLLZ);

NSPEC = size(gammax, 3);
% Find the element contains the source and its node
for ispec = 1 : NSPEC
    for j = 1:NGLLZ
        for i = 1:NGLLX
            if (ibool(i,j,ispec) == src_node)
                ispec_source = ispec;
                i_source = i;
                j_source = j;
                break;
            end
        end
    end
end
xi_source    = xigll(i_source);
gamma_source = zigll(j_source);

% Moment tensor solution
Mxx  =  Moment_Tensor(1,1);
Mzz  =  Moment_Tensor(3,3);
Mxz  =  Moment_Tensor(1,2);

[   hxis,    hpxis] = Lagrange_Any(   xi_source, NGLLX, xigll);
[hgammas, hpgammas] = Lagrange_Any(gamma_source, NGLLZ, zigll);

sourcearray = zeros(NDIM,NGLLX,NGLLZ);

if (strcmp(src_type, 'vector_source') == 1)
    % collocated force source
    for j = 1:NGLLZ
        for i = 1:NGLLX
            hlagrange = hxis(i) * hgammas(j);
            % source element is elastic, and only in P_SV case
            sourcearray(1,i,j) = -sin(force_angle*pi/180) * hlagrange;
            sourcearray(2,i,j) =  cos(force_angle*pi/180) * hlagrange;
        end
    end
elseif (strcmp(src_type, 'moment_tensor_source') == 1)
    dxis_dx    = 0.0;
    dxis_dz    = 0.0;
    dgammas_dx = 0.0;
    dgammas_dz = 0.0;
    for m = 1:NGLLZ
        for k = 1:NGLLX
            xixd    = xix(k,m,ispec_source);
            xizd    = xiz(k,m,ispec_source);
            gammaxd = gammax(k,m,ispec_source);
            gammazd = gammaz(k,m,ispec_source);
            
            hlagrange = hxis(k) * hgammas(m);
            
            dxis_dx    = dxis_dx + hlagrange * xixd;
            dxis_dz    = dxis_dz + hlagrange * xizd;
            dgammas_dx = dgammas_dx + hlagrange * gammaxd;
            dgammas_dz = dgammas_dz + hlagrange * gammazd;
        end
    end
    % calculate source array
    for m = 1:NGLLZ
        for k = 1:NGLLX
            dsrc_dx = (hpxis(k)*dxis_dx)*hgammas(m) + hpxis(k)*(hpgammas(m)*dgammas_dx);
            dsrc_dz = (hpxis(k)*dxis_dz)*hgammas(m) + hpxis(k)*(hpgammas(m)*dgammas_dz);
            sourcearray(1,k,m) = sourcearray(1,k,m)  + Mxx*dsrc_dx + Mxz*dsrc_dx;
            sourcearray(2,k,m) = sourcearray(2,k,m)  + Mxz*dsrc_dx + Mzz*dsrc_dz;
        end
    end
    
end

sourcearray = sourcearray / max(abs(sourcearray(:)));

end
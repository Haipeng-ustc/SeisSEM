function [Kx1, Kx2, Kx3, Kx4, Kx5, Kz1, Kz2, Kz3, Kz4, Kz5] = Assemble_Global_K(NGLOB, NGLLX, NGLLZ, jacobian, gammax, gammaz, xix, xiz, ibool, ...
                                                                                c11, c13, c33, c44, pml_dxx, pml_dzz)

%--------------------------------------
% Element wise local stiffness matrix:
%--------------------------------------
%% Stiffness matrix
%  K1  dphi_dx * dphi_dx
%  K2  dphi_dz * dphi_dz
%  K3  dphi_dx * dphi_dz
%  K4  dphi_dz * dphi_dx
%  K5  phi     * dphi_dx
%  K6  phi     * dphi_dz

    %%%%%%% Set up GLL points, weights and derivation matrices %%%%%%%%

    [xigll,wxgll,hprime_xx, hprimewgll_xx] = GetGLL(NGLLX);
    [zigll,wzgll,hprime_zz, hprimewgll_zz] = GetGLL(NGLLZ);

    NGLL = NGLLX;
    NSPEC = size(jacobian,3);
    % Stiffness matrix for temporary use
    K_ele_temp1 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    K_ele_temp2 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    K_ele_temp3 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    K_ele_temp4 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    K_ele_temp5 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    K_ele_temp6 = zeros(NGLLX,NGLLZ,NGLLX,NGLLZ);
    % Stiffness matrix for each element
    Ke1 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    Ke2 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    Ke3 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    Ke4 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    Ke5 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    Ke6 = zeros(NGLLX*NGLLZ,NGLLX*NGLLZ);
    % Stiffness Matrix assembly
    I  = zeros(length(ibool(:)),1);
    J  = zeros(length(ibool(:)),1);
    V1 = zeros(length(ibool(:)),1);
    V2 = zeros(length(ibool(:)),1);
    V3 = zeros(length(ibool(:)),1);
    V4 = zeros(length(ibool(:)),1);
    V5 = zeros(length(ibool(:)),1);
    V6 = zeros(length(ibool(:)),1);

    ct = 1;
    delta = eye(NGLL); % identity matrix: kronecker delta
    for ispec = 1:NSPEC
        K_ele_temp1(:) = 0.;
        K_ele_temp2(:) = 0.;
        K_ele_temp3(:) = 0.;
        K_ele_temp4(:) = 0.;
        K_ele_temp5(:) = 0.;
        K_ele_temp6(:) = 0.;

        for i = 1:NGLLZ
            for j = 1:NGLLX
                for k = 1:NGLLZ
                    for l = 1:NGLLX
                        term1 = 0.0; %  K1: phi_xx * phi_xx
                        term2 = 0.0; %  K2: phi_zz * phi_zz
                        term3 = 0.0; %  K3: phi_xx * phi_zz
                        term4 = 0.0; %  K4: phi_zz * phi_xx
                        term5 = 0.0; %  K5: phi    * phi_xx
                        term6 = 0.0; %  K6: phi    * phi_zz

                        for p = 1:NGLLX
                            for q = 1 : NGLLX
                                term1 = term1 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xix(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammax(p,q,ispec)) * ...
                                    (hprime_xx(k,p) * delta(l,q) * xix(p,q,ispec) + delta(k,p) * hprime_zz(l,q) * gammax(p,q,ispec)));
                                term2 = term2 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xiz(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammaz(p,q,ispec)) * ...
                                    (hprime_xx(k,p) * delta(l,q) * xiz(p,q,ispec) + delta(k,p) * hprime_zz(l,q) * gammaz(p,q,ispec)));
                                term3 = term3 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xix(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammax(p,q,ispec)) * ...
                                    (hprime_xx(k,p) * delta(l,q) * xiz(p,q,ispec) + delta(k,p) * hprime_zz(l,q) * gammaz(p,q,ispec)));
                                term4 = term4 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xiz(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammaz(p,q,ispec)) * ...
                                    (hprime_xx(k,p) * delta(l,q) * xix(p,q,ispec) + delta(k,p) * hprime_zz(l,q) * gammax(p,q,ispec)));
                                term5 = term5 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xix(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammax(p,q,ispec)) *  delta(k,p) *  delta(l,q));
                                term6 = term6 + wxgll(p) * wzgll(q) * jacobian(p,q,ispec) *  ...
                                    ((hprime_xx(i,p) * delta(j,q) * xiz(p,q,ispec) + delta(i,p) * hprime_zz(j,q) * gammaz(p,q,ispec)) *  delta(k,p) *  delta(l,q));
                            end
                        end
                        K_ele_temp1(i,j,k,l) = term1;
                        K_ele_temp2(i,j,k,l) = term2;
                        K_ele_temp3(i,j,k,l) = term3;
                        K_ele_temp4(i,j,k,l) = term4;
                        K_ele_temp5(i,j,k,l) = term5;
                        K_ele_temp6(i,j,k,l) = term6;

                    end
                end
            end
        end

        Ke1(:,:) = reshape(K_ele_temp1, NGLL*NGLL,NGLL*NGLL);
        Ke2(:,:) = reshape(K_ele_temp2, NGLL*NGLL,NGLL*NGLL);
        Ke3(:,:) = reshape(K_ele_temp3, NGLL*NGLL,NGLL*NGLL);
        Ke4(:,:) = reshape(K_ele_temp4, NGLL*NGLL,NGLL*NGLL);
        Ke5(:,:) = reshape(K_ele_temp5, NGLL*NGLL,NGLL*NGLL);
        Ke6(:,:) = reshape(K_ele_temp6, NGLL*NGLL,NGLL*NGLL);

        % Assembly: I,J,V are such that K(I,J) = V
        ig = ibool(:,:,ispec);
        ig = ig(:);
        for j = 1:length(ig)
            for i = 1:length(ig)
                I(ct)  = ig(i);
                J(ct)  = ig(j);
                V1(ct) = Ke1(i,j);
                V2(ct) = Ke2(i,j);
                V3(ct) = Ke3(i,j);
                V4(ct) = Ke4(i,j);
                V5(ct) = Ke5(i,j);
                V6(ct) = Ke6(i,j);

                ct = ct + 1;
            end
        end
    end
    % Faster assembly as SparseCSC format
    K1 = sparse(I,J,V1,NGLOB,NGLOB);
    K2 = sparse(I,J,V2,NGLOB,NGLOB);
    K3 = sparse(I,J,V3,NGLOB,NGLOB);
    K4 = sparse(I,J,V4,NGLOB,NGLOB);
    K5 = sparse(I,J,V5,NGLOB,NGLOB);
    K6 = sparse(I,J,V6,NGLOB,NGLOB);

    %% Stiffness matrix
    % Kx1 = c11 * dphi_dx * dphi_dx
    % Kx2 = c13 * dphi_dx * dphi_dz + c44 * dphi_dz * dphi_dx  ??
    % Kx3 = c44 * dphi_dz * dphi_dz
    % Kx4 = c11 * phi     * dphi_dx * pml_dxx
    % Kx5 = c44 * phi     * dphi_dz * pml_dzz
    % Kz1 = c44 * dphi_dx * dphi_dx
    % Kz2 = c44 * dphi_dx * dphi_dz + c13 * dphi_dz * dphi_dx  ??
    % Kz3 = c33 * dphi_dz * dphi_dz
    % Kz4 = c44 * phi     * dphi_dx * pml_dxx
    % Kz5 = c33 * phi     * dphi_dz * pml_dzz
    
    Kx1 = K1 .* c11;
    Kx2 = K3 .* c13 + K4 .* c44;
    Kx3 = K2 .* c44;
    Kx4 = K5 .* c11 .* pml_dxx;
    Kx5 = K6 .* c44 .* pml_dzz;

    Kz1 = K1 .* c44;
    Kz2 = K3 .* c44 + K4 .* c13;
    Kz3 = K2 .* c33;
    Kz4 = K5 .* c44 .* pml_dxx;
    Kz5 = K6 .* c33 .* pml_dzz;
    
end 
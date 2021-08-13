function [M, Cx, Cz, Cxx, Cxz, Czz] = Assemble_Global_M_and_C(NSPEC, NGLOB, NGLLX, NGLLZ, jacobian, ibool, rho, pml_dx, pml_dz)

    %% Mass matrix
    %       M   = rho * phi * phi
    %% Damping matrix
    %       Cx  = rho * phi * phi * pml_dx
    %       Cz  = rho * phi * phi * pml_dz
    %       Cxx = rho * phi * phi * pml_dx * pml_dx
    %       Cxz = rho * phi * phi * pml_dx * pml_dz
    %       Czz = rho * phi * phi * pml_dz * pml_dz
    
    [xigll,wxgll,hprime_xx, hprimewgll_xx] = GetGLL(NGLLX);
    [zigll,wzgll,hprime_zz, hprimewgll_zz] = GetGLL(NGLLZ);

    M    = zeros(NGLOB, 1);
    for ispec = 1:NSPEC
        for j = 1:NGLLZ
            for i = 1:NGLLX
                iglob = ibool(i,j,ispec);
                M(iglob)   =   M(iglob) + wxgll(i)* wzgll(j) * jacobian(i,j,ispec);
            end
        end
    end
    
    M   = M  .* rho;     
    Cx  = M  .* pml_dx;  
    Cz  = M  .* pml_dz;  
    Cxx = Cx .* pml_dx;  
    Cxz = Cx .* pml_dz;   
    Czz = Cz .* pml_dz;  
    
       
end

function  [coord, xix, xiz, gammax, gammaz, jacobian]= Assemble_Global_Jacobian(ibool,coorg,knods,NGLLX,NGLLZ)

if NGLLX ~=NGLLZ
    error('NGLLX is not the same with NGLLZ');
end 


[xigll,wxgll,hprime_xx, hprimewgll_xx] = GetGLL(NGLLX);
[zigll,wzgll,hprime_zz, hprimewgll_zz] = GetGLL(NGLLZ);

% These parameters are the same with the main function.
ngnod = size(knods,1);
NSPEC = size(ibool,3);
NDIM  = size(coorg,1);
NGLOB = max(ibool(:));

coord = zeros(NDIM,NGLOB); % set the coordinates of the points of the global grid
xix   = zeros(NGLLX, NGLLZ, NSPEC);
xiz   = zeros(NGLLX, NGLLZ, NSPEC);
gammax   = zeros(NGLLX, NGLLZ, NSPEC);
gammaz   = zeros(NGLLX, NGLLZ, NSPEC);
jacobian = zeros(NGLLX, NGLLZ, NSPEC);

for ispec = 1:NSPEC
    for j = 1:NGLLZ
        for i = 1:NGLLX
            
            xi    = xigll(i);
            gamma = zigll(j);
            
            % recompute jacobian for any (xi,gamma) point, not necessarily a GLL point
            
            % create the 2D shape functions and then the Jacobian
            [shape2D, dershape2D] = Shape_Functions(xi,gamma,ngnod,NDIM);
            
            % compute coordinates and jacobian matrix
            x = 0.0;
            z = 0.0;
            
            xxi = 0.0;
            zxi = 0.0;
            xgamma = 0.0;
            zgamma = 0.0;

            for ia=1:ngnod

                nnum = knods(ia,ispec);

                xelm = coorg(1,nnum);
                zelm = coorg(2,nnum);

                x = x + shape2D(ia)*xelm;
                z = z + shape2D(ia)*zelm;

                xxi = xxi + dershape2D(1,ia)*xelm;
                zxi = zxi + dershape2D(1,ia)*zelm;
                xgamma = xgamma + dershape2D(2,ia)*xelm;
                zgamma = zgamma + dershape2D(2,ia)*zelm;

            end
            jacobianl = xxi*zgamma - xgamma*zxi;
      
            if (jacobianl <= 0.d0)
                % If the Jacobian is negative, it means that there is an error in the mesh
                % print the coordinates of the mesh points of this element
                fprintf('ispec = %d\n', ispec);
                fprintf('ngnod = %d\n', ngnod);
                for ia=1:ngnod
                    nnum = knods(ia,ispec);
                    xelm = coorg(1,nnum);
                    zelm = coorg(2,nnum);
                    fprintf('node %d, at numeber %d, (x,y) = %f, %f\n', ia,nnum,xelm,zelm);
                end
                error('error: negative 2D Jacobian found');
            end

            % invert the relation
            xixl    =   zgamma / jacobianl;
            gammaxl = - zxi    / jacobianl;
            xizl    = - xgamma / jacobianl;
            gammazl =   xxi    / jacobianl;

            % assign to global array
            coord(1,ibool(i,j,ispec)) = x;
            coord(2,ibool(i,j,ispec)) = z;
            
            xix(i,j,ispec) = xixl;
            xiz(i,j,ispec) = xizl;
            gammax(i,j,ispec) = gammaxl;
            gammaz(i,j,ispec) = gammazl;
            jacobian(i,j,ispec) = jacobianl;
            
        end
    end
end

end
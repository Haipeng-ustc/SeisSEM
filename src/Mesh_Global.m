function [ibool] = Mesh_Global(nelem_x, nelem_z, NGLOB, NGLL, NSPEC)
% clear;clc
% 
% NGLL = 5;
% NGLLX = NGLL;
% NGLLZ = NGLL;
% 
% nelem_x = 6;
% nelem_z = 8;
% 
% NSPEC = nelem_x * nelem_z;                                 % number of spectral elements of the mesh
% NGLOB = (nelem_x*(NGLLX-1) + 1) * (nelem_z*(NGLLZ-1) + 1); % number of unique grid points of the mesh
% NPGEO = (nelem_x+1)*(nelem_z+1);   


ibool = zeros(NGLL,NGLL,NSPEC);	        % local to global index mapping
igL  = reshape([1:NGLL*(NGLL-1)],NGLL-1,NGLL);
igB  = reshape([1:NGLL*(NGLL-1)],NGLL,NGLL-1);
igLB = reshape([1:(NGLL-1)*(NGLL-1)],NGLL-1,NGLL-1);
e    = 0;
last_iglob = 0;

% use the nelem_z in ex direction 
% use the nelem_x in ez direction, it is wired!  
for ex=1:nelem_z
    for ez=1:nelem_x
        e = e+1;
        % Take care of redundant nodes at element edges :
        if e==1               % first element: bottom-left
            ig = reshape([1:NGLL*NGLL],NGLL,NGLL);
        else
            if ex==1 	      % elements on first (left) column
                ig(:,     1) = ibool(:,NGLL,e-1); 		% bottom edge
                ig(:,2:NGLL) = last_iglob + igB;        % the rest
            elseif ez==1      % elements on first (bottom) row
                ig(1,     :) = ibool(NGLL,:,e-nelem_x); % left edge
                ig(2:NGLL,:) = last_iglob + igL; 		% the rest
            else              % other elements
                ig(1,:) = ibool(NGLL,:,e-nelem_x); 		% left edge
                ig(:,1) = ibool(:,NGLL,e-1); 		    % bottom edge
                ig(2:NGLL,2:NGLL) = last_iglob + igLB;
            end
        end
        ibool(:,:,e) = ig;
        last_iglob = ig(NGLL,NGLL);
    end
end

%% change here
ibool = permute(ibool, [2  1  3]);


end






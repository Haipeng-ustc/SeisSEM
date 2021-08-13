function [coorg, knods, ibool] = Read_Gmsh(meshpath, NGLLX, NGLLZ)

    if NGLLX ~=5 && NGLLZ~=5
        error('Only 5 points element is supported, while other cases are under development');
    end 
    % open the generated Gmsh
    fp    = fopen(meshpath,'r');
    for i = 1:4
        temp = fgetl(fp);
    end
    node_num = fscanf(fp, '%d\n',1);
    coorg  = zeros(node_num,2);
    for  i = 1:node_num
        value = fscanf(fp, '%d %f %f %f\n',4);
        if value(1) ~= i
            error('Wrong node number!');
        end
        coorg(i,1) = value(2);
        coorg(i,2) = value(3);
    end
    for i = 1:2
        temp = fgetl(fp);
    end
    element_num  = fscanf(fp, '%d\n',1);
    ibool = zeros(element_num,NGLLX*NGLLZ);
    iele = 0;
    for  i = 1:element_num
        value = fgetl(fp);
        value = str2double(strsplit(value,' '));
        if length(value) == NGLLX*NGLLZ + 5     % 5 addtional values
            iele = iele+1;
            ibool(iele,1:NGLLX*NGLLZ) = value(end-NGLLX*NGLLZ+1:end); 
        end
    end
    element_num  = iele;
    ibool = ibool(1:element_num,:);
    
    if min(ibool(:)) ~=1 || max(ibool(:)) ~=node_num
        error('Element node numerbing is wrong!');
    end
    
    % Geometry node numbering 
    knods = fliplr(ibool(:,1:4)); % 4-node control point, clockwise numerbing 
       
    % Adjust global GLL nodes numerbing
    index_change_25_nodes = [4, 14, 15, 16, 1, 13, 20, 24, 17, 5, 12, 23, 25, 21, 6, 11, 19, 22, 18, 7, 3, 10, 9, 8, 2];
    for i = 1 :  element_num
        ibool(i,:) = ibool(i,index_change_25_nodes);
    end
    ibool = reshape(ibool', [NGLLX, NGLLZ, element_num]);
    ibool = permute(ibool, [2,1,3]);
    
    % Take the transpose
    knods = knods';
    coorg = coorg';
        
end


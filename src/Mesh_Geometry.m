function [node_xz, element_node] = Mesh_Geometry(xmin, xmax, zmin, zmax, Nx, Nz, type)

%     I  ELEMENT_TYPE   Definition
%     -  ------------   ----------
%     1  T3             3 node linear triangle;
%     2  T6             6 node quadratic triangle;
%     3  T10            10 node cubic triangle.
%     4  Q4             4 node linear Lagrange/serendipity quadrilateral;
%     5  Q9             9 node quadratic Lagrange quadrilateral;
%     6  Q16            16 node cubic Lagrange quadrilateral;

% Geometrical points that describe the mesh:
%        node_xz: coordinates of all the corners of the mesh elements
%   element_node: numbering of the four corners of each mesh element



if ~strcmp(type, 'q4' )
    error('Set type as q4, other mesh is under construction ...');
end 

if (strcmp(type, 't3' ))
    
    %
    % T3  Grid:
    %
    %     9---10---11---12
    %     |\ 8 |\10 |\12 |
    %     | \  | \  | \  |
    %     |  \ |  \ |  \ |
    %     |  7\|  9\| 11\|
    %     5----6----7----8
    %     |\ 2 |\ 4 |\ 6 |
    %     | \  | \  | \  |
    %     |  \ |  \ |  \ |
    %     |  1\|  3\|  5\|
    %     1----2----3----4
    %
    Node_num    = (Nx + 1) * (Nz + 1);
    
    Element_num = 2 * Nx * Nz;
    
    element = 1;
    for j = 1 : Nz
        for i = 1 : Nx
            
            sw = i     + ( j - 1 ) * ( Nx + 1 );
            se = i + 1 + ( j - 1 ) * ( Nx + 1 );
            nw = i     +   j       * ( Nx + 1 );
            ne = i + 1 +   j       * ( Nx + 1 );
            
            element_node(element, 1) = sw;
            element_node(element, 2) = se;
            element_node(element, 3) = nw;
            element = element + 1;
            
            element_node(element, 1) = ne;
            element_node(element, 2) = nw;
            element_node(element, 3) = se;
            element = element + 1;
            
        end
    end
    
    for j = 1 : Nz + 1
        for i = 1 : Nx + 1
            
            node_xz(i + (j-1) * (Nx + 1), 1 ) = ((Nx - i + 1) * xmin + ( i - 1 )*xmax) / (1.0 * Nx);
            node_xz(i + (j-1) * (Nx + 1), 2 ) = ((Nz - j + 1) * zmin + ( j - 1 )*zmax) / (1.0 * Nz);
            
        end
    end
    
    
elseif(strcmp(type, 't6' ))
    
    %
    %
    % T6  Grid:
    %
    %    29-30-31-32-33-34-35
    %     |\ 8  |\10  |\12  |
    %     | \   | \   | \   |
    %    22 23 24 25 26 27 28
    %     |   \ |   \ |   \ |
    %     |  7 \|  9 \| 11 \|
    %    15-16-17-18-19-20-21
    %     |\ 2  |\ 4  |\ 6  |
    %     | \   | \   | \   |
    %     8  9 10 11 12 13 14
    %     |   \ |   \ |   \ |
    %     |  1 \|  3 \|  5 \|
    %     1--2--3--4--5--6--7
    %
    
    Node_num = (2 * Nx + 1) * (2 * Nz + 1);
    
    Element_num = 2 * Nx * Nz;
    
    element = 1;
    
    for j = 1 : Nz
        for i = 1 : Nx
            
            sw = 2 * ( j - 1 )  * ( 2 * Nx + 1 ) + 2 * ( i - 1 ) + 1;
            w  = sw +               2 * Nx + 1;
            nw = sw +         2 * ( 2 * Nx + 1 );
            
            s  = sw + 1;
            c  = sw + 1 +           2 * Nx + 1;
            n  = sw + 1 +     2 * ( 2 * Nx + 1 );
            
            se = sw + 2;
            e  = sw + 2 +           2 * Nx + 1;
            ne = sw + 2 +     2 * ( 2 * Nx + 1 );
            
            element_node(element, 1) = sw;
            element_node(element, 2) = se;
            element_node(element, 3) = nw;
            element_node(element, 4) = s;
            element_node(element, 5) = c;
            element_node(element, 6) = w;
            element = element + 1;
            
            element_node(element, 1) = ne;
            element_node(element, 2) = nw;
            element_node(element, 3) = se;
            element_node(element, 4) = n;
            element_node(element, 5) = c;
            element_node(element, 6) = e;
            element = element + 1;
            
        end
    end
    
    for j = 1 : 2 * Nz + 1
        for i = 1 : 2 * Nx + 1
            
            node_xz( i + (j - 1) * (2 * Nx + 1), 1) = (i - 1) * (xmax - xmin) / (2.0 * Nx);
            node_xz( i + (j - 1) * (2 * Nx + 1), 2) = (j - 1) * (zmax - zmin) / (2.0 * Nz);
            
        end
    end
    
    
    
elseif(strcmp(type, 't10' ))
    
    
    %
    % T10 Grid:
    %
    %    43-44-45-46-47-48-49
    %     |\     6 |\     8 |
    %     | \      | \      |
    %    36 37 38 39 40 41 42
    %     |   \    |   \    |
    %     |    \   |    \   |
    %    29 30 31 32 33 34 35
    %     |      \ |      \ |
    %     | 5     \| 7     \|
    %    22-23-24-25-26-27-28
    %     |\     2 |\     4 |
    %     | \      | \      |
    %    15 16 17 18 19 20 21
    %     |   \    |   \    |
    %     |    \   |    \   |
    %     8  9 10 11 12 13 14
    %     |      \ |      \ |
    %     | 1     \| 3     \|
    %     1--2--3--4--5--6--7
    %
    
    Node_num = (3 * Nx + 1) * (3 * Nz + 1);
    
    Element_num = 2 * Nx * Nz;
    
    element = 1;
    
    for j = 1 : Nz
        for i = 1 : Nx
            
            base = (j - 1) * 3 * (3 * Nx + 1) + 3 * i - 2;
            
            element_node(element, 1) = base;
            element_node(element, 2) = base + 1;
            element_node(element, 3) = base + 2;
            element_node(element, 4) = base + 3;
            element_node(element, 5) = base + (3 * Nx + 1) + 2;
            element_node(element, 6) = base + 2 * (3 * Nx + 1) + 1;
            element_node(element, 7) = base + 3 * (3 * Nx + 1);
            element_node(element, 8) = base + 2 * (3 * Nx + 1);
            element_node(element, 9) = base + (3 * Nx + 1);
            element_node(element, 10) = base + (3 * Nx + 1) + 1;
            element = element + 1;
            
            element_node(element, 1) = base + 3 * (3 * Nx + 1) + 3;
            element_node(element, 2) = base + 3 * (3 * Nx + 1) + 2;
            element_node(element, 3) = base + 3 * (3 * Nx + 1) + 1;
            element_node(element, 4) = base + 3 * (3 * Nx + 1);
            element_node(element, 5) = base + 2 * (3 * Nx + 1) + 1;
            element_node(element, 6) = base + (3 * Nx + 1) + 2;
            element_node(element, 7) = base + 3;
            element_node(element, 8) = base + (3 * Nx + 1) + 3;
            element_node(element, 9) = base + 2 * (3 * Nx + 1) + 3;
            element_node(element, 10) = base + 2 * (3 * Nx + 1) + 2;
            element = element + 1;
            
        end
    end
    
    for j = 1 : 3 * Nz + 1
        for i = 1 : 3 * Nx + 1
            
            node_xz( i + (j - 1) * (3 * Nx + 1), 1) = (i - 1) * (xmax - xmin) / (3.0 * Nx);
            node_xz( i + (j - 1) * (3 * Nx + 1), 2) = (j - 1) * (zmax - zmin) / (3.0 * Nz);
            
        end
    end
    
    
elseif(strcmp(type, 'q4' ))
    
    %
    %   Q4  Grid:
    %
    %     9---10---11---12
    %     |    |    |    |
    %     |    |    |    |
    %     |  4 |  5 |  6 |
    %     |    |    |    |
    %     5----6----7----8
    %     |    |    |    |
    %     |    |    |    |
    %     |  1 |  2 |  3 |
    %     |    |    |    |
    %     1----2----3----4
    %
    
    Node_num = (Nx + 1) * (Nz + 1);
    
    Element_num =  Nx * Nz;
    
    element = 1;
    
    for  j = 1 : Nz
        for i = 1 : Nx
            
            sw = i     + ( j - 1 ) * ( Nx + 1 );
            se = i + 1 + ( j - 1 ) * ( Nx + 1 );
            nw = i     +   j       * ( Nx + 1 );
            ne = i + 1 +   j       * ( Nx + 1 );
            
            element_node(element, 1) = sw;
            element_node(element, 2) = se;
            element_node(element, 3) = ne;
            element_node(element, 4) = nw;
            
            element = element + 1;
            
        end
    end
    
    for j = 1 : Nz + 1
        for i = 1 : Nx + 1
            
            node_xz(i + ( j - 1 ) * (Nx + 1),1 ) = ((Nx - i + 1) * xmin + (i - 1) * xmax) / (1.0 * Nx);
            node_xz(i + ( j - 1 ) * (Nx + 1),2 ) = ((Nz - j + 1) * zmin + (j - 1) * zmax) / (1.0 * Nz);
            
        end
    end
    
    
elseif(strcmp(type, 'q9' ))
    
    
    %
    %
    %    Q9    Grid:
    %
    %    29---30---31---32---33---34---35
    %     |    .    |    .    |    .    |
    %     |    .    |    .    |    .    |
    %    22 . 23 . 24 . 25 . 26 . 27 . 28
    %     |    .    |    .    |    .    |
    %     | 4  .    | 5  .    | 6  .    |
    %    15---16---17---18---19---20---21
    %     |    .    |    .    |    .    |
    %     |    .    |    .    |    .    |
    %     8 .  9 . 10 . 11 . 12 . 13 . 14
    %     |    .    |    .    |    .    |
    %     | 1  .    | 2  .    | 3  .    |
    %     1----2----3----4----5----6----7
    %
    %
    
    Node_num = (2 * Nx + 1) * (2 * Nz + 1);
    
    Element_num =  Nx * Nz;
    
    element = 1;
    
    for j = 1 : Nz
        
        for i = 1 : Nx
            
            sw = 2 * ( j - 1 )  * ( 2 * Nx + 1 ) + 2 * ( i - 1 ) + 1;
            w  = sw +               2 * Nx + 1;
            nw = sw +         2 * ( 2 * Nx + 1 );
            
            s  = sw + 1;
            c  = sw + 1 +               2 * Nx + 1;
            n  = sw + 1 +         2 * ( 2 * Nx + 1 );
            
            se = sw + 2;
            e  = sw + 2 +               2 * Nx + 1;
            ne = sw + 2 +         2 * ( 2 * Nx + 1 );
            
            element_node(element, 1) = sw;
            element_node(element, 2) = se;
            element_node(element, 3) = ne;
            element_node(element, 4) = nw;
            element_node(element, 5) = s;
            element_node(element, 6) = e;
            element_node(element, 7) = n;
            element_node(element, 8) = w;
            element_node(element, 9) = c;
            
            element = element + 1;
            
        end
    end
    
    for j = 1 : 2 * Nz + 1
        for i = 1 : 2 * Nx + 1
            
            node_xz(i + (j - 1 ) * (2 * Nx + 1), 1 ) = ( i - 1 ) * (xmax - xmin) / (2.0 * Nx);
            node_xz(i + (j - 1 ) * (2 * Nx + 1), 2 ) = ( j - 1 ) * (zmax - zmin) / (2.0 * Nz);
            
        end
    end
    
elseif(strcmp(type, 'q16' ))
    
    %
    %   q16  Grid:
    %
    %    43-44-45-46-47-48-49
    %     |        |        |
    %     |        |        |
    %    36 37 38 39 40 41 42
    %     |        |        |
    %     |        |        |
    %    29 30 31 32 33 34 35
    %     |        |        |
    %     | 3      | 4      |
    %    22-23-24-25-26-27-28
    %     |        |        |
    %     |        |        |
    %    15 16 17 18 19 20 21
    %     |        |        |
    %     |        |        |
    %     8  9 10 11 12 13 14
    %     |        |        |
    %     | 1      | 2      |
    %     1--2--3--4--5--6--7
    %
    
    Node_num = (3 * Nx + 1) * (3 * Nz + 1);
    
    Element_num =  Nx * Nz;
    
    
    element = 1;
    
    for j = 1 : Nz
        
        for i = 1 : Nx
            
            base = (j - 1) * 3 * (3 * Nx + 1) + 3 * i - 2;
            
            element_node(element, 1) = base;
            element_node(element, 2) = base + 1;
            element_node(element, 3) = base + 2;
            element_node(element, 4) = base + 3;
            element_node(element, 5) = base + (3 * Nx + 1);
            element_node(element, 6) = base + (3 * Nx + 1) + 1;
            element_node(element, 7) = base + (3 * Nx + 1) + 2;
            element_node(element, 8) = base + (3 * Nx + 1) + 3;
            element_node(element, 9) = base + 2 * (3 * Nx + 1);
            element_node(element, 10) = base + 2 * (3 * Nx + 1) + 1;
            element_node(element, 11) = base + 2 * (3 * Nx + 1) + 2;
            element_node(element, 12) = base + 2 * (3 * Nx + 1) + 3;
            element_node(element, 13) = base + 3 * (3 * Nx + 1);
            element_node(element, 14) = base + 3 * (3 * Nx + 1) + 1;
            element_node(element, 15) = base + 3 * (3 * Nx + 1) + 2;
            element_node(element, 16) = base + 3 * (3 * Nx + 1) + 3;
            
            element = element + 1;
            
        end
    end
    for  j = 1 : 3 * Nz + 1
        for i = 1 : 3 * Nx + 1
            
            node_xz(i + (j - 1 ) * (3 * Nx + 1), 1 ) = ( i - 1 ) * (xmax - xmin) / (3.0 * Nx);
            node_xz(i + (j - 1 ) * (3 * Nx + 1), 2 ) = ( j - 1 ) * (zmax - zmin) / (3.0 * Nz);
            
        end
    end
    
else
    
    Node_num = -1;
    Element_num = -1;
    node_xz    = -1;
    element_node = -1;
    fprintf('Generate_Mesh - Fatal error!\n');
    fprintf('Wrong element type: %s\n', type);
    return;
    
end

node_xz = node_xz';
element_node = element_node';


end






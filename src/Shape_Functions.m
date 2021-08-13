function [shape2D, dershape2D] = Shape_Functions(xi,gamma,ngnod,NDIM)
%=======================================================================
%
%                               4 . . . . 7 . . . . 3
%                               .                   .
%                               .         gamma     .
%                               .                   .
%                               8         9  xi     6
%                               .                   .
%                               .                   .
%                               .                   .
%                               1 . . . . 5 . . . . 2
%
%                           Local coordinate system : s,t
%
%=======================================================================

shape2D    = zeros(1,ngnod);
dershape2D = zeros(NDIM,ngnod);

% very small value
TINYVAL = 1.e-9;

%---- set up the shape functions and their local derivatives
s  = xi;
t  = gamma;

%----    4-node element
if (ngnod == 4)
    sp = s + 1.0;
    sm = s - 1.0;
    tp = t + 1.0;
    tm = t - 1.0;
    
    %----  corner nodes
    shape2D(1) =   0.25 * sm * tm;
    shape2D(2) = - 0.25 * sp * tm;
    shape2D(3) =   0.25 * sp * tp;
    shape2D(4) = - 0.25 * sm * tp;
    
    dershape2D(1,1) =   0.25 * tm;
    dershape2D(1,2) = - 0.25 * tm;
    dershape2D(1,3) =   0.25 * tp;
    dershape2D(1,4) = - 0.25 * tp;
    
    dershape2D(2,1) =   0.25 * sm;
    dershape2D(2,2) = - 0.25 * sp;
    dershape2D(2,3) =   0.25 * sp;
    dershape2D(2,4) = - 0.25 * sm;
    
%----    9-node element
elseif (ngnod == 9)
    
    sp = s + 1;
    sm = s - 1;
    tp = t + 1;
    tm = t - 1;
    s2 = s * 2;
    t2 = t * 2;
    ss = s * s;
    tt = t * t;
    st = s * t;
    
    %----  corner nodes
    shape2D(1) = 0.25 * sm * st * tm;
    shape2D(2) = 0.25 * sp * st * tm;
    shape2D(3) = 0.25 * sp * st * tp;
    shape2D(4) = 0.25 * sm * st * tp;
    
    dershape2D(1,1) = 0.25 * tm * t * (s2 - 1.0);
    dershape2D(1,2) = 0.25 * tm * t * (s2 + 1.0);
    dershape2D(1,3) = 0.25 * tp * t * (s2 + 1.0);
    dershape2D(1,4) = 0.25 * tp * t * (s2 - 1.0);
    
    dershape2D(2,1) = 0.25 * sm * s * (t2 - 1.0);
    dershape2D(2,2) = 0.25 * sp * s * (t2 - 1.0);
    dershape2D(2,3) = 0.25 * sp * s * (t2 + 1.0);
    dershape2D(2,4) = 0.25 * sm * s * (t2 + 1.0);
    
    %----  midside nodes
    shape2D(5) = 0.5 * tm * t * (1.0 - ss);
    shape2D(6) = 0.5 * sp * s * (1.0 - tt);
    shape2D(7) = 0.5 * tp * t * (1.0 - ss);
    shape2D(8) = 0.5 * sm * s * (1.0 - tt);
    
    dershape2D(1,5) = -1.0 * st * tm;
    dershape2D(1,6) =  0.5 * (1.0 - tt) * (s2 + 1.0);
    dershape2D(1,7) = -1.0 * st * tp;
    dershape2D(1,8) =  0.5 * (1.0 - tt) * (s2 - 1.0);
    
    dershape2D(2,5) =  0.5 * (1.0 - ss) * (t2 - 1.0);
    dershape2D(2,6) = -1.0 * st * sp;
    dershape2D(2,7) =  0.5 * (1.0 - ss) * (t2 + 1.0);
    dershape2D(2,8) = -1.0 * st * sm;
    
    %----  center node
    shape2D(9) = (1.0 - ss) * (1.0 - tt);
   
    dershape2D(1,9) = -1.0 * s2 * (1.0 - tt);
    dershape2D(2,9) = -1.0 * t2 * (1.0 - ss);
    
else
    error ('Error: wrong number of control nodes');
end

%--- check the shape functions and their derivatives
% sum of shape functions should be one
if (abs(sum(shape2D)-1.0) > TINYVAL)
    error('error shape functions');
end
% sum of derivatives of shape functions should be zero
if (abs(sum(dershape2D(1,:))) > TINYVAL)
    error( 'error deriv xi shape functions')
end
if (abs(sum(dershape2D(2,:))) > TINYVAL)
    error('error deriv gamma shape functions')
end


end
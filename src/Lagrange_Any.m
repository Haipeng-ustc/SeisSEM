function [h,hprime] = Lagrange_Any(xi,NGLL,xigll)
% subroutine to compute the Lagrange interpolants based upon the interpolation points
% and their first derivatives at any point xi in [-1,1]

h      = zeros(NGLL, 1);
hprime = zeros(NGLL, 1);

for dgr = 1:NGLL
    prod1 = 1.0;
    prod2 = 1.0;
    % lagrangian interpolants
    x0 = xigll(dgr);
    for i = 1:NGLL
        if (i ~= dgr)
            x = xigll(i);
            prod1 = prod1*(xi-x);
            prod2 = prod2*(x0-x);
        end
    end
    
    % takes inverse to avoid additional divisions (multiplications are cheaper than divisions)
    prod2_inv = 1.0/prod2;
    h(dgr)    = prod1 * prod2_inv;
    
    % first derivatives
    sum = 0.0;
    for i = 1: NGLL
        if (i ~= dgr)
            prod3 = 1.0;
            for j = 1:NGLL
                if (j ~= dgr && j ~= i)
                    prod3 = prod3*(xi-xigll(j));
                end
                sum = sum + prod3;
            end
        end
        hprime(dgr) = sum * prod2_inv;
    end
end
end 
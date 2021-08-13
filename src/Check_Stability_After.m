function Check_Stability_After(displ, NGLOB, it, deltat, NSTEP, STABILITY_THRESHOLD)

    Usolidnorm = -1.0;
    for iglob = 1:NGLOB
        current_value = sqrt(displ(1,iglob)^2 + displ(2,iglob)^2);
        if (current_value > Usolidnorm)
            Usolidnorm = current_value;
        end
    end
    fprintf('Time step %6d out of %d\n',it,NSTEP);
    % compute current time
    time = (it-1)*deltat;
    % check stability of the code, exit if unstable
    if (Usolidnorm > STABILITY_THRESHOLD || Usolidnorm < 0)
        fprintf('Max norm of displacement is = %f\n',Usolidnorm);
        error('The SEM code became unstable and blew up.');
    end
end
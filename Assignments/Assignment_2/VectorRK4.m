function [t_rk4, Y_rk4] = VectorRK4(ode_function, t0, Y0, NStep, h, a, b, c, d, e)
    % Pre-allocation 
    t_rk4 = zeros(NStep, 1);
    Y_rk4 = zeros(NStep, length(Y0));

    % to store results of y1 and y2 in a singe array
    Y_rk4(1, :) = Y0; 
    
    % initial condition
    t_rk4(1) = t0; 

    for k = 1:NStep
        s1 = ode_function(t_rk4(k), Y_rk4(k, :), a, b, c, d, e);
        s2 = ode_function(t_rk4(k) + h/2, Y_rk4(k, :) + h/2 .* s1', a, b, c, d, e);
        s3 = ode_function(t_rk4(k) + h/2, Y_rk4(k, :) + h/2 .* s2', a, b, c, d, e);
        s4 = ode_function(t_rk4(k) + h, Y_rk4(k, :) + h .* s3', a, b, c, d, e);
        Y_rk4(k + 1, :) = Y_rk4(k, :) + h * (s1/6 + s2/3 + s3/3 + s4/6)';
        t_rk4(k + 1) = t_rk4(k) + h;
    end
end

function [t, Y] = VectorRK4(ode_system, tspan, Y0, NStep, h)
    % Initialize variables
    t = linspace(tspan(1), tspan(2), NStep+1);
    Y = zeros(NStep+1, numel(Y0));
    Y(1, :) = Y0;

    % Perform Runge-Kutta integration
    for k = 1:NStep
        % Compute the slopes
        s1 = ode_system(t(k), Y(k, :)')';
        s2 = ode_system(t(k) + h/2, (Y(k, :)' + h/2 * s1)')';
        s3 = ode_system(t(k) + h/2, (Y(k, :)' + h/2 * s2)')';
        s4 = ode_system(t(k) + h, (Y(k, :)' + h * s3)')';

        % Update the solution
        Y(k+1, :) = Y(k, :) + h * (s1/6 + s2/3 + s3/3 + s4/6);
    end
end

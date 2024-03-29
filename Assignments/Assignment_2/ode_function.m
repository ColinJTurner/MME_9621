function dydt = ode_function(t, Y, a, b, c, d, e)
    % Extract variables
    x1 = Y(1);
    v1 = Y(2);
    x2 = Y(3);
    v2 = Y(4);
    x3 = Y(5);
    v3 = Y(6);

    % Compute derivatives
    dx1dt = v1;
    dv1dt = -a .* x1 + b .* (x2 - x1);
    dx2dt = v2;
    dv2dt = c .* (x1 - x2) + d .* (x3 - x2);
    dx3dt = v3;
    dv3dt = e .* (x2 - x3);

    % Return the derivatives
    dydt = [dx1dt; dv1dt; dx2dt; dv2dt; dx3dt; dv3dt];
end
% function dydt = circuit_odes(t, y, V, R, L1, L2, M, C1, C2)
%     % Extract variables
%     i1 = y(1);
%     vc1 = y(2);
%     i2 = y(3);
%     vc2 = y(4);
% 
%     % Define system of ODEs
%     dydt = zeros(4,1);
%     dydt(1) = (V(t) - vc1 - R*i1 - M/L1*i2) / L1;
%     dydt(2) = i1 / C1;
%     dydt(3) = (vc1 - vc2 + M/L2*i1 - R*i2) / L2;
%     dydt(4) = i2 / C2;
% end

function dydt = circuit_odes(t, y, V, R, L1, L2, M, C1, C2)
    % Extract variables
    i1 = y(1);
    vc1 = y(2);
    i2 = y(3);
    vc2 = y(4);

    % Define system of ODEs
    dydt = zeros(4,1);
    dydt(1) = (V(t) - vc1 - R*i1 + M*dydt(3)) / L1;
    dydt(2) = i1 / C1;
    dydt(3) = (V(t) - vc2 - R*i2 + M*dydt(1)) / L2;
    dydt(4) = i2 / C2;
end

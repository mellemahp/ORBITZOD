function [state_acc] = Full_Nonlinear_Dynamics(C, t, state_m, noise) 
    % Generates accelarations based on full nonlinear EOM
    state_acc(1) = state_m(2);
    state_acc(2) = -C.mu * state_m(1) / sqrt(state_m(1)^2 + state_m(3)^2)^3 + noise(1);
    state_acc(3) = state_m(4);
    state_acc(4) = -C.mu * state_m(3) / sqrt(state_m(1)^2 + state_m(3)^2)^3 + noise(2);
    state_acc = state_acc';
end

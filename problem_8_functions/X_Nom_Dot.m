function [x_nom_dot_out] = X_Nom_Dot(C, t)
    % Calculates the nominal x velocity at a given time
    x_nom_dot_out = -C.r0 * C.n * sin(C.n * t); 
end

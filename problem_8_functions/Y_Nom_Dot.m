function [y_nom_dot_out] = Y_Nom_Dot(C, t)
    % Calculates the nominal y velocity at a given time
    y_nom_dot_out = C.r0 * C.n * cos(C.n * t);
end 


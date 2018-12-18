function [y_nom_out] = Y_Nom(C, t) 
    % Calculates the nominal y position at a given time
    y_nom_out = C.r0 * sin(C.n * t);
end 

function [x_nom_out] = X_Nom(C, t)
    % Calculates the nominal x position at a given time
    x_nom_out = C.r0 * cos(C.n * t); 
end


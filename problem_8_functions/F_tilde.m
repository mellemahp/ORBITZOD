function [f_tilde_out] = F_tilde(C, t)
    %==========================================================================
    % [f_tilde_out] = F_tilde(C, time)
    %
    % Consructs the linearized dynamics matrix (linearized about the
    % nominal trajectory) evaluated at a given time
    % time
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to evaluate F_tilde           sec
    %
    % OUTPUT:       zeros(size(xp_0,1), size(msrs_nom,1))'; 
    %    
    %  f_tilde_out      - taylor expanded dynamics matrix for S/C        (4x4)
    %                                     
    % Coupling:
    %   
    %  None 
    %
    %==========================================================================
    
    f_tilde_out = eye(4) + C.delta_t * [0, 1, 0, 0;
                                        C.mu * (3 * (cos(C.n * t))^2 - 1) / (C.r0^3), 0, ...
                                        3 * C.mu * sin(2 * C.n * t)/(2 * C.r0^3), 0;
                                        0, 0, 0, 1;
                                        3 * C.mu * sin(2 * C.n * t) / (2 * C.r0^3), 0, ...
                                        C.mu * ( 3 * (sin(C.n * t))^2 - 1) / (C.r0^3), 0];
end
 
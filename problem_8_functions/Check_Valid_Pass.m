

  function bool = Check_Valid_Pass(C, t, i)
    %==========================================================================
    % [bool] = check_valid_pass(C, t, i)
    %
    % Determines whether a spacecraft is visible within a given stations
    % viewing cone
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  t                   - time at which to initialize h matrix        sec
    %  i                   - number of station to generate h matrix for  int [1-12] 
    %
    % OUTPUT:       
    %    
    %  bool               - boolean as to whether S/C is visible to stn  (boolean)
    %                                     
    % Coupling:
    %
    %  phi_i              - finds elevation angle of spacecraft as viewed
    %                       from station
    %  theta_i            - finds angle of station in inertial frame
    %
    %==========================================================================
    u = [X_Nom(C, t) - X_i(C, t, i), Y_Nom(C,t) - Y_i(C, t, i), 0];
    v = [X_i(C, t, i), Y_i(C, t, i), 0];
    
    ThetaBetween = atan2(norm(cross(u,v)),dot(u,v));
    
    if -pi/2 < ThetaBetween && ThetaBetween < pi/2
        bool = true;
    else 
        bool = false;
    end 
 end


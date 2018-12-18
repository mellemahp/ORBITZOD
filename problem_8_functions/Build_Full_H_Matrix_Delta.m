 function [h_matrix] = Build_Full_H_Matrix_Delta(C, time)
    %==========================================================================
    % [h_matrix] = Build_Full_H_Matrix(C, time)
    %
    % Consructs the perturbation observation matrix for all of the stations at a 
    % given time
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %
    % OUTPUT:       
    %    
    %  h_matrix      - matrix of gps satellite orbit parameters         (nx25)
    %                                     
    % Coupling:
    %
    % Check Valid Pass => Uses this to check that sat is within valid
    % ground station viewing cone
    %
    % H_i => Uses this to construct the observation for each individual ground 
    %        station  
    %
    %==========================================================================
    h_matrix = [];
    for i = 1:12
        if Check_Valid_Pass(C, time, i)
            h_matrix = [h_matrix; H_i(C, time, i)];
        else
            h_matrix = [h_matrix; NaN(3, 4)];
        end
    end   
 end 

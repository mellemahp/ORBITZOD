function [h_matrix] = Get_Msrs_True(C, sc_state, time, R) 
    h_matrix = [];
    for i = 1:12
        if Check_Valid_Pass_Truth(C, sc_state, time, i)
            h_new = True_Msrs(C, sc_state, time, i) + mvnrnd([0, 0, 0], R, 1)';
            h_matrix = [h_matrix; h_new];
        else
            h_matrix = [h_matrix; NaN(3, 1)];
        end
    end   
end 


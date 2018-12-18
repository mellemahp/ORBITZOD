function Plot_Measurements(times, msrs)
    % Plots the Range, Rangerate, and Phi measurements for a given dataset
    figure()

    x = 1;
    for i = 1:12
        subplot(3, 1, 1)
        hold on 
        scatter(times, msrs(x,1:length(times)))
        ylabel('Range (km)')
        subplot(3, 1, 2) 
        hold on 
        scatter(times, msrs(x + 1,1:length(times)))
        ylabel('Range Rate (km/s)')
        subplot(3, 1, 3)
        hold on
        scatter(times, msrs(x + 2, 1:length(times)))
        ylabel('Phi (rad)')
        xlabel('Time (sec)')
        x = x + 3; 
    end
end

function msrs_corrected = Make_Data_Useful(data)
    msrs_corrected = NaN(36, length(data)); 
    index_list = 1:3:36;
    for i = 1:length(data)
        if ~isempty(data{i})
            for j = 1:length(data{i}(1,:))
                stn_dat = data{i}(:,j); 
                if ~isnan(stn_dat(4))
                    stn_id = stn_dat(4); 
                    stn_idx = index_list(stn_id);
                    msrs_corrected(stn_idx:stn_idx+2,i) = stn_dat(1:3);
                end
            end
        end
    end
end

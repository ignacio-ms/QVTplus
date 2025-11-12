function [Locs]=get_SampleLocs(data_struct,Labels)
    %made up get_SampleLocs
    % Initialize Locs array
    Locs = -1 * ones(1, length(Labels));
    
    % Loop through each label in Labels and find the corresponding location
    for i = 1:length(Labels)
        loc = Labels{i, 3}; % Assuming locs is the 3rd column in Labels
        n_loc = str2num(loc);
        n_lab = str2num(Labels{i,2});

        if sum(size(n_loc)) > 2
            subset = find(data_struct.branchList(:,4) == n_loc(1));
            idx = subset(n_loc(2));

        elseif sum(size(n_loc)) == 2
            subset = find(data_struct.branchList(:,4) == n_lab);
            idx = subset(n_loc);
        end
        % Find the matching index in data_struct.branchList based on the location
        
        % If found, store the index in Locs
        if ~isempty(idx)
            Locs(i) = idx;
        end
    end

end
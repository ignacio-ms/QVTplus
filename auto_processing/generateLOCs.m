function [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT)
    % Generate the LOCs structure based on the data_struct and correspondenceDict
    % Handles ICA, BA cutoff logic, and main vessel LOC identification

    LOCs = struct(); % Initialize the LOCs structure

    %% Step 1: Process ICA and BA vessels
    ica_slice = 0;
    LICA = correspondenceDict.LICA;
    RICA = correspondenceDict.RICA;
    BA = correspondenceDict.BASI;

    % Extract branch information for each vessel
    info_LICA = find_LOCs('extractBranchInfo', data_struct, LICA);
    info_RICA = find_LOCs('extractBranchInfo', data_struct, RICA);
    info_BA = find_LOCs('extractBranchInfo', data_struct, BA);

    % Round Z values
    info_LICA(:, 3) = round(info_LICA(:, 3));
    info_RICA(:, 3) = round(info_RICA(:, 3));
    info_BA(:, 3) = round(info_BA(:, 3));

    % Find unique Z values
    unique_z_LICA = unique(info_LICA(:, 3));
    unique_z_RICA = unique(info_RICA(:, 3));
    unique_z_BA = unique(info_BA(:, 3));

    % Find the largest Z slice with one point for each vessel
    max_single_z = find_LOCs('findMaxSingleZ', unique_z_LICA, unique_z_RICA, unique_z_BA, ...
                             info_LICA, info_RICA, info_BA);

    if max_single_z > -inf
        ica_slice = max_single_z;

        % Extract candidate locations around the slice (±1 Z tolerance)
        % This allows us to select from multiple candidates based on circularity
        z_tolerance = 5;  % Allow ±1 Z slice for candidate selection
        LICA_candidates = info_LICA(abs(info_LICA(:, 3) - max_single_z) <= z_tolerance, :);
        RICA_candidates = info_RICA(abs(info_RICA(:, 3) - max_single_z) <= z_tolerance, :);
        BA_candidates = info_BA(abs(info_BA(:, 3) - max_single_z) <= z_tolerance, :);
        
        % Select best LOC based on circularity
        LICA_LOC = find_LOCs('selectLOCByCircularity', LICA_candidates, data_struct, info_LICA, 4, 2);
        RICA_LOC = find_LOCs('selectLOCByCircularity', RICA_candidates, data_struct, info_RICA, 4, 2);
        BA_LOC = find_LOCs('selectLOCByCircularity', BA_candidates, data_struct, info_BA, 4, 2);

        % Store the LOCs for ICA and BA
        if ~isempty(LICA_LOC) && size(LICA_LOC, 2) >= 5
            LOCs.LICA = [LICA_LOC(1, 4), LICA_LOC(1, 5)];
        end
        if ~isempty(RICA_LOC) && size(RICA_LOC, 2) >= 5
            LOCs.RICA = [RICA_LOC(1, 4), RICA_LOC(1, 5)];
        end
        if ~isempty(BA_LOC) && size(BA_LOC, 2) >= 5
            LOCs.BASI = [BA_LOC(1, 4), BA_LOC(1, 5)];
        end
    end
    %% Step 2: Handle venous system
    vesselLabels = fieldnames(correspondenceDict);
    for i = 1:numel(vesselLabels)
        keyName = vesselLabels{i};
        if strcmp(keyName, 'SSSV')
            SSSV = correspondenceDict.SSSV;
            SSSV_LOC = find_LOCs('extractSSSV',SSSV,data_struct);
            try
                LOCs.SSSV = [SSSV_LOC(1, 4), SSSV_LOC(1, 5)];
            catch
                continue
            end
        end
        if strcmp(keyName, 'LTSV')
            LTSV = correspondenceDict.LTSV;
            LTSV_LOC = find_LOCs('extractLTSV',LTSV,data_struct);
            try
                LOCs.LTSV = [LTSV_LOC(1, 4), LTSV_LOC(1, 5)];
            catch
                continue
            end
        end
        if strcmp(keyName, 'RTSV')
            RTSV = correspondenceDict.RTSV;
            RTSV_LOC = find_LOCs('extractRTSV',RTSV,data_struct);
            try
                LOCs.RTSV = [RTSV_LOC(1, 4), RTSV_LOC(1, 5)];
            catch
                continue
            end
        end
        if strcmp(keyName, 'STRV')
            STRV = correspondenceDict.STRV;
            % Check if SSSV was already assigned and get its segment ID
            sssv_segment_id = [];
            if isfield(LOCs, 'SSSV') && ~isempty(LOCs.SSSV)
                sssv_segment_id = LOCs.SSSV(1);
            end
            STRV_LOC = find_LOCs('extractSTRV',STRV,data_struct, sssv_segment_id);
            try
                if ~isempty(STRV_LOC)
                    LOCs.STRV = [STRV_LOC(1, 4), STRV_LOC(1, 5)];
                end
            catch
                continue
            end
        end
    end

    if isfield(LOCs, 'SSSV') && isfield(LOCs, 'STRV')
        if LOCs.SSSV(1) == LOCs.STRV(1)
            LOCs = rmfield(LOCs, 'STRV');
        end
    end
    %TODO: revisa que L/R esta bien.
    LOCs = resolveLongVenousSegment(LOCs, data_struct, 'RTSV');
    LOCs = resolveLongVenousSegment(LOCs, data_struct, 'LTSV');


    %% Step 3: Process other vessels
    for i = 1:numel(vesselLabels)
        keyName = vesselLabels{i};

        if strcmp(keyName, 'LICA') || strcmp(keyName, 'RICA') || strcmp(keyName, 'BASI')
            % Already processed ICA and BA
            continue;
        elseif ismember(keyName, {'LPCA', 'RPCA', 'LMCA', 'RMCA', 'RACA', 'LACA'})
            % Process main vessels
            LOCs.(keyName) = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT);
        elseif strcmp(keyName, 'COMM')
            % Special case for COMM vessels
            LOCs.(keyName) = unique(correspondenceDict.COMM);
        end
    end

    %% Step 4: Handle PCA and secondary PCA logic
    if isfield(correspondenceDict, 'RPC2') && isfield(LOCs, 'RPCA')
        if ismember(LOCs.('RPCA')(1), correspondenceDict.('RPC2'))
            correspondenceDict.('RPCA') = [correspondenceDict.('RPCA'); correspondenceDict.('RPC2')];
            correspondenceDict = rmfield(correspondenceDict, 'RPC2');
        else
            correspondenceDict = rmfield(correspondenceDict, 'RPC2');
        end
    end
    
    if isfield(correspondenceDict, 'LPC2') && isfield(LOCs, 'LPCA')
        if ismember(LOCs.('LPCA')(1), correspondenceDict.('LPC2'))
            correspondenceDict.('LPCA') = [correspondenceDict.('LPCA'); correspondenceDict.('LPC2')];
            correspondenceDict = rmfield(correspondenceDict, 'LPC2');
        else
            correspondenceDict = rmfield(correspondenceDict, 'LPC2');
        end
    end
end


function LOCs = resolveLongVenousSegment(LOCs, data_struct, fieldName)
    if isfield(LOCs, 'SSSV') && isfield(LOCs, fieldName)
        if LOCs.SSSV(1) == LOCs.(fieldName)(1)
            shared_seg_id = LOCs.SSSV(1);
            branchList = data_struct.branchList;
            all_points = branchList(branchList(:, 4) == shared_seg_id, :);
            n_points = size(all_points, 1);

            step = floor(n_points / 6);
            parts = cell(6,1);
            for i = 1:6
                idx_start = (i-1)*step + 1;
                if i < 6
                    idx_end = i*step;
                else
                    idx_end = n_points;
                end
                parts{i} = all_points(idx_start:idx_end, :);
            end

            % Most vertical (FIRST 3)
            best_vert_score = -Inf;
            best_vert_idx = NaN;
            for i = 1:3
                p = parts{i};
                score = std(p(:,3));
                if score > best_vert_score
                    best_vert_score = score;
                    best_vert_idx = i;
                end
            end
            vert_part = parts{best_vert_idx};
            mid_idx = floor(size(vert_part,1)/2);
            SSSV_point = vert_part(mid_idx, :);
            [~, SSSV_row_idx] = min(sum(abs(branchList - SSSV_point),2));
            LOCs.SSSV = [shared_seg_id, SSSV_row_idx];

            % Most horizontal (LAST 3)
            best_horiz_score = Inf;
            best_horiz_idx = NaN;
            for i = 4:6
                p = parts{i};
                score = std(p(:,3));  % low std(Z) → horizontal
                if score < best_horiz_score
                    best_horiz_score = score;
                    best_horiz_idx = i;
                end
            end
            horiz_part = parts{best_horiz_idx};
            mid_idx = floor(size(horiz_part,1)/2);
            target_point = horiz_part(mid_idx, :);
            [~, row_idx] = min(sum(abs(branchList - target_point),2));
            LOCs.(fieldName) = [shared_seg_id, row_idx];
        end
    end
end


function [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT)
    % Generate the LOCs structure based on the data_struct and correspondenceDict
    % Handles ICA, BA cutoff logic, and main vessel LOC identification

    LOCs = struct(); % Initialize the LOCs structure

    %% Step 1: Process ICA and BA vessels (skip if any key missing, e.g. centerline_test with different segments)
    fprintf('[generateLOCs] Step 1 (ICA/BA): correspondenceDict has LICA=%d, RICA=%d, BASI=%d\n', ...
        isfield(correspondenceDict, 'LICA'), isfield(correspondenceDict, 'RICA'), isfield(correspondenceDict, 'BASI'));
    ica_slice = 0;
    if isfield(correspondenceDict, 'LICA'), LICA = correspondenceDict.LICA; else, LICA = []; end
    if isfield(correspondenceDict, 'RICA'), RICA = correspondenceDict.RICA; else, RICA = []; end
    if isfield(correspondenceDict, 'BASI'), BA = correspondenceDict.BASI; else, BA = []; end

    % Extract branch information for each vessel
    info_LICA = find_LOCs('extractBranchInfo', data_struct, LICA);
    info_RICA = find_LOCs('extractBranchInfo', data_struct, RICA);
    info_BA = find_LOCs('extractBranchInfo', data_struct, BA);
    fprintf('[generateLOCs] Branch info points: LICA=%d, RICA=%d, BA=%d\n', ...
        size(info_LICA, 1), size(info_RICA, 1), size(info_BA, 1));

    % Round Z values (only when non-empty to avoid subscript errors)
    if ~isempty(info_LICA), info_LICA(:, 3) = round(info_LICA(:, 3)); end
    if ~isempty(info_RICA), info_RICA(:, 3) = round(info_RICA(:, 3)); end
    if ~isempty(info_BA), info_BA(:, 3) = round(info_BA(:, 3)); end

    % Find unique Z values
    if ~isempty(info_LICA), unique_z_LICA = unique(info_LICA(:, 3)); else, unique_z_LICA = []; end
    if ~isempty(info_RICA), unique_z_RICA = unique(info_RICA(:, 3)); else, unique_z_RICA = []; end
    if ~isempty(info_BA), unique_z_BA = unique(info_BA(:, 3)); else, unique_z_BA = []; end

    % Find the largest Z slice with one point for each vessel
    max_single_z = find_LOCs('findMaxSingleZ', unique_z_LICA, unique_z_RICA, unique_z_BA, ...
                             info_LICA, info_RICA, info_BA);
    fprintf('[generateLOCs] max_single_z (ICA/BA slice) = %g (valid if > -inf)\n', max_single_z);

    if max_single_z > -inf
        ica_slice = max_single_z;

        % Extract candidate locations around the slice (±1 Z tolerance)
        % This allows us to select from multiple candidates based on circularity
        z_tolerance = 5;  % Allow ±1 Z slice for candidate selection
        LICA_candidates = info_LICA(abs(info_LICA(:, 3) - max_single_z) >= z_tolerance, :);
        RICA_candidates = info_RICA(abs(info_RICA(:, 3) - max_single_z) >= z_tolerance, :);
        BA_candidates = info_BA(abs(info_BA(:, 3) - max_single_z) >= z_tolerance, :);

        % Select best LOC based on circularity
        fprintf('[generateLOCs] Selecting LICA LOC based on circularity\n');
        LICA_LOC = find_LOCs('selectICAsBASILOC', LICA_candidates, data_struct);
        fprintf('[generateLOCs] Selecting RICA LOC based on circularity\n');
        RICA_LOC = find_LOCs('selectICAsBASILOC', RICA_candidates, data_struct);
        fprintf('[generateLOCs] Selecting BA LOC based on circularity\n');
        BA_LOC = find_LOCs('selectICAsBASILOC', BA_candidates, data_struct);
        fprintf('[generateLOCs] ICA/BA candidates: LICA=%d, RICA=%d, BA=%d; LOC selected: LICA=%d, RICA=%d, BA=%d\n', ...
            size(LICA_candidates, 1), size(RICA_candidates, 1), size(BA_candidates, 1), ...
            ~isempty(LICA_LOC) && size(LICA_LOC, 2) >= 5, ~isempty(RICA_LOC) && size(RICA_LOC, 2) >= 5, ~isempty(BA_LOC) && size(BA_LOC, 2) >= 5);

        % Store the LOCs for ICA and BA
        if ~isempty(LICA_LOC) && size(LICA_LOC, 2) >= 5
            LOCs.LICA = [LICA_LOC(1, 4), LICA_LOC(1, 5)];
        else
            fprintf('[generateLOCs] LICA LOC not stored (empty or <5 cols)\n');
        end
        if ~isempty(RICA_LOC) && size(RICA_LOC, 2) >= 5
            LOCs.RICA = [RICA_LOC(1, 4), RICA_LOC(1, 5)];
        else
            fprintf('[generateLOCs] RICA LOC not stored (empty or <5 cols)\n');
        end
        if ~isempty(BA_LOC) && size(BA_LOC, 2) >= 5
            LOCs.BASI = [BA_LOC(1, 4), BA_LOC(1, 5)];
        else
            fprintf('[generateLOCs] BASI LOC not stored (empty or <5 cols)\n');
        end
    else
        fprintf('[generateLOCs] ICA/BA skipped: no common Z slice (max_single_z = -inf). Check that LICA, RICA, BA all have points.\n');
    end
    %% Step 2: Handle venous system
    vesselLabels = fieldnames(correspondenceDict);
    for i = 1:numel(vesselLabels)
        keyName = vesselLabels{i};
        if strcmp(keyName, 'SSSV')
            SSSV = correspondenceDict.SSSV;
            try
                SSSV_LOC = find_LOCs('extractSSSV',SSSV,data_struct);
                LOCs.SSSV = [SSSV_LOC(1, 4), SSSV_LOC(1, 5)];
                fprintf('[generateLOCs] Venous SSSV: LOC stored [seg %d, pt %d]\n', LOCs.SSSV(1), LOCs.SSSV(2));
            catch ME
                fprintf('[generateLOCs] Venous SSSV: skipped (%s)\n', ME.message);
            end
        end
        if strcmp(keyName, 'LTSV')
            LTSV = correspondenceDict.LTSV;
            try
                LTSV_LOC = find_LOCs('extractLTSV',LTSV,data_struct);
                LOCs.LTSV = [LTSV_LOC(1, 4), LTSV_LOC(1, 5)];
                fprintf('[generateLOCs] Venous LTSV: LOC stored [seg %d, pt %d]\n', LOCs.LTSV(1), LOCs.LTSV(2));
            catch ME
                fprintf('[generateLOCs] Venous LTSV: skipped (%s)\n', ME.message);
            end
        end
        if strcmp(keyName, 'RTSV')
            RTSV = correspondenceDict.RTSV;
            try
                RTSV_LOC = find_LOCs('extractRTSV',RTSV,data_struct);
                LOCs.RTSV = [RTSV_LOC(1, 4), RTSV_LOC(1, 5)];
                fprintf('[generateLOCs] Venous RTSV: LOC stored [seg %d, pt %d]\n', LOCs.RTSV(1), LOCs.RTSV(2));
            catch ME
                fprintf('[generateLOCs] Venous RTSV: skipped (%s)\n', ME.message);
            end
        end
        if strcmp(keyName, 'STRV')
            STRV = correspondenceDict.STRV;
            sssv_segment_id = [];
            if isfield(LOCs, 'SSSV') && ~isempty(LOCs.SSSV)
                sssv_segment_id = LOCs.SSSV(1);
            end
            try
                STRV_LOC = find_LOCs('extractSTRV',STRV,data_struct, sssv_segment_id);
                if ~isempty(STRV_LOC)
                    LOCs.STRV = [STRV_LOC(1, 4), STRV_LOC(1, 5)];
                    fprintf('[generateLOCs] Venous STRV: LOC stored [seg %d, pt %d]\n', LOCs.STRV(1), LOCs.STRV(2));
                else
                    fprintf('[generateLOCs] Venous STRV: skipped (empty STRV_LOC)\n');
                end
            catch ME
                fprintf('[generateLOCs] Venous STRV: skipped (%s)\n', ME.message);
            end
        end
    end

    % When SSSV and another venous share the same segment, resolve to distinct points.
    %TODO: revisa que L/R esta bien.
    LOCs = resolveLongVenousSegment(LOCs, data_struct, 'RTSV');
    LOCs = resolveLongVenousSegment(LOCs, data_struct, 'LTSV');
    LOCs = resolveSSSVSTRV(LOCs, data_struct);
    LOCs = validateSSSVSTRVSwap(LOCs, data_struct);

    fprintf('[generateLOCs] Venous LOCs after resolving long venous segments:\n');
    if isfield(LOCs, 'SSSV')
        fprintf('[generateLOCs]   SSSV: [seg %d, pt %d]\n', LOCs.SSSV(1), LOCs.SSSV(2));
    end
    if isfield(LOCs, 'RTSV')
        fprintf('[generateLOCs]   RTSV: [seg %d, pt %d]\n', LOCs.RTSV(1), LOCs.RTSV(2));
    end
    if isfield(LOCs, 'LTSV')
        fprintf('[generateLOCs]   LTSV: [seg %d, pt %d]\n', LOCs.LTSV(1), LOCs.LTSV(2));
    end
    if isfield(LOCs, 'STRV')
        fprintf('[generateLOCs]   STRV: [seg %d, pt %d]\n', LOCs.STRV(1), LOCs.STRV(2));
    end

    %% Step 3: Process other vessels
    fprintf('[generateLOCs] Step 3 (main vessels): processing %s\n', strjoin(vesselLabels, ', '));
    for i = 1:numel(vesselLabels)
        keyName = vesselLabels{i};

        if strcmp(keyName, 'LICA') || strcmp(keyName, 'RICA') || strcmp(keyName, 'BASI')
            continue;
        elseif ismember(keyName, {'LPCA', 'RPCA', 'LMCA', 'RMCA', 'RACA', 'LACA'})
            locResult = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT);
            if ~isempty(locResult) && numel(locResult) >= 2 && isfinite(locResult(1)) && isfinite(locResult(2))
                LOCs.(keyName) = locResult(1:2);
                fprintf('[generateLOCs]   %s: LOC stored [seg %d, pt %d]\n', keyName, locResult(1), locResult(2));
            else
                fprintf('[generateLOCs]   %s: LOC missing (processMainVessels returned empty/invalid)\n', keyName);
            end
        elseif strcmp(keyName, 'RCOMM') || strcmp(keyName, 'LCOMM')
            % Process RCOMM/LCOMM vessels (split from COMM based on RAS orientation)
            % Use similar logic to main vessels: find best LOC in the segment
            % Handle missing vessels like venous vessels (use try-catch and check validity)
            if isfield(correspondenceDict, keyName) && ~isempty(correspondenceDict.(keyName))
                try
                    locResult = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT);
                    % Check if result is valid (not NaN and has 2 elements)
                    if ~isempty(locResult) && numel(locResult) == 2 && ...
                       ~isnan(locResult(1)) && ~isnan(locResult(2)) && ...
                       isfinite(locResult(1)) && isfinite(locResult(2))
                        LOCs.(keyName) = locResult;
                        disp('RCOMM/LCOMM vessel found and added to LOCs');
                    end
                    % If invalid, don't add to LOCs (vessel doesn't exist or couldn't be found)
                catch
                    % If extraction fails, skip this vessel (like venous vessels)
                    disp('RCOMM/LCOMM vessel not found and skipped');
                    continue
                end
            end
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

    % Debug: summary of assigned LOCs
    locKeys = fieldnames(LOCs);
    expectedLoc = {'LICA', 'RICA', 'BASI', 'LMCA', 'RMCA', 'LPCA', 'RPCA', 'LACA', 'RACA', 'RCOMM', 'LCOMM', 'SSSV', 'LTSV', 'RTSV', 'STRV'};
    missingLoc = setdiff(expectedLoc, locKeys);
    fprintf('[generateLOCs] Summary: %d LOCs assigned. Missing: %s\n', numel(locKeys), strjoin(missingLoc, ', '));
end


function LOCs = resolveLongVenousSegment(LOCs, data_struct, fieldName)
    if isfield(LOCs, 'SSSV') && isfield(LOCs, fieldName)
        if LOCs.SSSV(1) == LOCs.(fieldName)(1)
            shared_seg_id = LOCs.SSSV(1);
            branchList = data_struct.branchList;
            all_points = branchList(branchList(:, 4) == shared_seg_id, :);
            parts = splitIntoParts(all_points, 6);

            % Most vertical (FIRST 3)
            best_vert_score = -Inf;
            best_vert_idx = NaN;
            for i = 1:3
                score = std(parts{i}(:, 3));
                if score > best_vert_score
                    best_vert_score = score;
                    best_vert_idx = i;
                end
            end
            SSSV_row = pickMaskedMidpoint(parts{best_vert_idx}, branchList, data_struct);
            LOCs.SSSV = [shared_seg_id, SSSV_row];

            % Most horizontal (LAST 3)
            best_horiz_score = Inf;
            best_horiz_idx = NaN;
            for i = 4:6
                score = std(parts{i}(:, 3));
                if score < best_horiz_score
                    best_horiz_score = score;
                    best_horiz_idx = i;
                end
            end
            other_row = pickMaskedMidpoint(parts{best_horiz_idx}, branchList, data_struct);
            LOCs.(fieldName) = [shared_seg_id, other_row];
        end
    end
end


function LOCs = resolveSSSVSTRV(LOCs, data_struct)
    % When SSSV and STRV share the same segment, resolve to distinct points using
    % principal direction (SVD): STRV is more diagonal/vertical (~45-90°, aligned with
    % [0,1,1]), SSSV is more horizontal (runs along AP at the top of the brain).
    % This is orientation-independent unlike a simple Y split.
    if ~isfield(LOCs, 'SSSV') || ~isfield(LOCs, 'STRV')
        return;
    end
    if LOCs.SSSV(1) ~= LOCs.STRV(1)
        return;
    end
    shared_seg_id = LOCs.SSSV(1);
    branchList = data_struct.branchList;
    all_points = branchList(branchList(:, 4) == shared_seg_id, :);
    parts = splitIntoParts(all_points, 6);

    % STRV expected direction: diagonal in Y-Z plane (~45°), same as extractSTRV
    strv_expected = [0; 1; 1];
    strv_expected = strv_expected / norm(strv_expected);

    % Score each part by alignment with STRV direction (higher = more STRV-like)
    nParts = numel(parts);
    alignment = zeros(nParts, 1);
    for i = 1:nParts
        pts = parts{i}(:, 1:3);
        if size(pts, 1) < 3
            alignment(i) = 0;
            continue;
        end
        centered = pts - mean(pts);
        [~, ~, V] = svd(centered, 'econ');
        direction = V(:, 1);
        alignment(i) = abs(dot(direction, strv_expected));
    end

    % STRV = most aligned with diagonal direction; SSSV = least aligned (most horizontal)
    [~, strv_idx] = max(alignment);
    [~, sssv_idx] = min(alignment);

    % Avoid assigning both to the same part
    if strv_idx == sssv_idx
        sorted_align = sort(alignment);
        sssv_idx = find(alignment == sorted_align(1), 1);
        strv_idx = find(alignment == sorted_align(end), 1);
    end

    SSSV_clIdx = pickMaskedMidpoint(parts{sssv_idx}, branchList, data_struct);
    STRV_clIdx = pickMaskedMidpoint(parts{strv_idx}, branchList, data_struct);
    LOCs.SSSV = [shared_seg_id, SSSV_clIdx];
    LOCs.STRV = [shared_seg_id, STRV_clIdx];
end


function LOCs = validateSSSVSTRVSwap(LOCs, data_struct)
    % Check if SSSV and STRV LOCs are swapped by comparing the local centerline
    % direction at each LOC against expected anatomy:
    %   STRV: diagonal/vertical (~45-90°), aligned with [0,1,1]
    %   SSSV: more horizontal (runs along AP at top of brain), low alignment with [0,1,1]
    % If the SSSV LOC sits in a STRV-like direction and vice versa, swap them.
    if ~isfield(LOCs, 'SSSV') || ~isfield(LOCs, 'STRV')
        return;
    end
    branchList = data_struct.branchList;
    strv_expected = [0; 1; 1];
    strv_expected = strv_expected / norm(strv_expected);

    sssv_align = localDirectionAlignment(LOCs.SSSV, branchList, strv_expected);
    strv_align = localDirectionAlignment(LOCs.STRV, branchList, strv_expected);

    % STRV should have HIGHER alignment with [0,1,1] than SSSV
    if sssv_align > strv_align
        fprintf('[validateSSSVSTRVSwap] Swapping SSSV and STRV (SSSV alignment=%.3f > STRV alignment=%.3f)\n', ...
            sssv_align, strv_align);
        tmp = LOCs.SSSV;
        LOCs.SSSV = LOCs.STRV;
        LOCs.STRV = tmp;
    end
end


function align = localDirectionAlignment(loc, branchList, expected_dir)
    % Compute alignment of the local centerline direction at a LOC point with an expected direction.
    % Uses a window of ±5 neighboring points along the same segment to estimate direction via SVD.
    seg = loc(1);
    clIdx = loc(2);
    segRows = find(branchList(:, 4) == seg);
    rowMatch = find(branchList(:, 4) == seg & branchList(:, 5) == clIdx, 1);
    if isempty(rowMatch)
        align = 0;
        return;
    end
    pos = find(segRows == rowMatch);
    if isempty(pos)
        align = 0;
        return;
    end
    window = 5;
    lo = max(1, pos - window);
    hi = min(numel(segRows), pos + window);
    neighborRows = segRows(lo:hi);
    pts = branchList(neighborRows, 1:3);
    if size(pts, 1) < 3
        align = 0;
        return;
    end
    centered = pts - mean(pts);
    [~, ~, V] = svd(centered, 'econ');
    direction = V(:, 1);
    align = abs(dot(direction, expected_dir));
end


function parts = splitIntoParts(all_points, nParts)
    n_points = size(all_points, 1);
    step = floor(n_points / nParts);
    parts = cell(nParts, 1);
    for i = 1:nParts
        idx_start = (i - 1) * step + 1;
        if i < nParts
            idx_end = i * step;
        else
            idx_end = n_points;
        end
        parts{i} = all_points(idx_start:idx_end, :);
    end
end


function row_idx = pickMaskedMidpoint(part_points, branchList, data_struct)
    % From a segment part, pick the point closest to the midpoint that is inside the
    % 4D flow segmentation mask. Falls back to geometric midpoint if none are in-mask.
    segVol = data_struct.segment;
    dims = size(segVol);
    n = size(part_points, 1);
    inMask = false(n, 1);
    for i = 1:n
        iy = round(part_points(i, 1));
        ix = round(part_points(i, 2));
        iz = round(part_points(i, 3));
        if iy >= 1 && iy <= dims(1) && ix >= 1 && ix <= dims(2) && iz >= 1 && iz <= dims(3)
            inMask(i) = segVol(iy, ix, iz) > 0;
        end
    end

    if any(inMask)
        masked_pts = part_points(inMask, :);
        mid_idx = floor(size(masked_pts, 1) / 2);
        mid_idx = max(1, mid_idx);
        chosen_point = masked_pts(mid_idx, :);
    else
        mid_idx = max(1, floor(n / 2));
        chosen_point = part_points(mid_idx, :);
    end
    [~, row_idx] = min(sum(abs(branchList - chosen_point), 2));
end


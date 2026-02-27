function varargout = find_LOCs(funcName, varargin)
    % Dispatcher function
    switch funcName
        case 'extractBranchInfo'
            varargout{1} = extractBranchInfo(varargin{:});
        case 'findMaxSingleZ'
            varargout{1} = findMaxSingleZ(varargin{:});
        case 'extractLocation'
            varargout{1} = extractLocation(varargin{:});
        case 'ensureMinZOffset'
            varargout{1} = ensureMinZOffset(varargin{:});
        case 'extractSSSV'
            varargout{1} = extractSSSV(varargin{:});
        case 'extractLTSV'
            varargout{1} = extractLTSV(varargin{:});
        case 'extractRTSV'
            varargout{1} = extractRTSV(varargin{:});
        case 'extractSTRV'
            varargout{1} = extractSTRV(varargin{:});
        case 'getCircularity'
            varargout{1} = getCircularity(varargin{:});
        case 'selectICAsBASILOC'
            varargout{1} = selectICAsBASILOC(varargin{:});
        otherwise
            error('Unknown function name: %s', funcName);
    end
end

function branchInfo = extractBranchInfo(data_struct, branchLabels)
    branchInfo = [];
    for i = 1:length(branchLabels)
        branchInfo = [branchInfo; data_struct.branchList(data_struct.branchList(:, 4) == branchLabels(i), :)];
    end
end

function maxZ = findMaxSingleZ(z_LICA, z_RICA, z_BA, info_LICA, info_RICA, info_BA)
    maxZ = -inf;
    for z = max([z_LICA; z_RICA; z_BA]):-1:min([z_LICA; z_RICA; z_BA])
        count_LICA = sum(int32(info_LICA(:, 3)) == z);
        count_RICA = sum(int32(info_RICA(:, 3)) == z);
        count_BA = sum(int32(info_BA(:, 3)) == z);

        if count_LICA >= 1 && count_RICA >= 1 && count_BA >= 1
            maxZ = z;
            break;
        end
    end
end

function loc = extractLocation(info, zValue, zColumn)
    loc = info(find(info(:, zColumn) == zValue, 1), :);
end

function loc = ensureMinZOffset(info, loc, minOffset, topOffset)
    % Ensure minimum Z offset is at least minOffset
    if loc(1, 5) < minOffset
        loc(1, 5) = minOffset;
    end

    % Calculate max Z value for the same branch
    maxZ = max(info(info(:, 4) == loc(4), 5));

    % Adjust if the difference is less than minOffset
    if (maxZ - loc(1, 5)) < topOffset
        loc(1, 5) = maxZ - topOffset;
    end
end

function loc = extractSSSV(segment_ids, data_struct)
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);

    % Only consider points inside the 4D flow segmentation mask
    inMask = isInsideSegmentMask(data_struct, branchList);
    branchList = branchList(inMask, :);
    if isempty(branchList), loc = []; return; end

    x = branchList(:,1);
    y = branchList(:,2);
    z = branchList(:,3);
    segment_id = branchList(:,4);

    dims = size(data_struct.segment);
    minX = min(x); maxX = max(x); 
    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);
    min_thrX = minX + 0.33 * (maxX - minX);
    max_thrX = maxX - 0.33 * (maxX - minX);
    posterior_thresh = minY + 0.4 * (maxY - minY);
    superior_thresh = minZ + 0.33 * (maxZ - minZ);
    posterior_mask = y <= posterior_thresh;
    superior_mask  = z >= superior_thresh;

    ss_mask = posterior_mask & superior_mask;
    subset = branchList(ss_mask, :);

    unique_segments = unique(subset(:,4));
    num_segments = length(unique_segments);

    if num_segments == 0
        loc = [];
        return;
    elseif num_segments <= 2
        counts = arrayfun(@(s) sum(subset(:,4) == s), unique_segments);
        [~, idx] = max(counts);
        ss_segment_id = unique_segments(idx);
    else
        lengths = zeros(num_segments,1);
        for i = 1:num_segments
            seg_id = unique_segments(i);
            lengths(i) = sum(subset(:,4) == seg_id);
        end
        [~, idx] = max(lengths);
        ss_segment_id = unique_segments(idx);
    end

    seg_points = branchList(segment_id == ss_segment_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);  % Return full row
    end

function loc = extractLTSV(segment_ids, data_struct)
    loc = extractLateralTSV(segment_ids, data_struct, 'left');
end

function loc = extractRTSV(segment_ids, data_struct)
    loc = extractLateralTSV(segment_ids, data_struct, 'right');
end

function loc = extractLateralTSV(segment_ids, data_struct, side)
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);

    % Only consider points inside the 4D flow segmentation mask
    inMask = isInsideSegmentMask(data_struct, branchList);
    branchList = branchList(inMask, :);
    if isempty(branchList), loc = []; return; end

    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); val = branchList(:,4);

    minX = min(x); maxX = max(x);
    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);

    % x_left_thresh  = minX + 0.33 * (maxX - minX);
    % x_right_thresh = maxX - 0.33 * (maxX - minX);
    % posterior_thresh = minY + 0.33 * (maxY - minY);
    % z_inferior_thresh = minZ + 0.7 * (maxZ - minZ);
    x_left_thresh  = minX + 0.40 * (maxX - minX);
    x_right_thresh = maxX - 0.40 * (maxX - minX);
    posterior_thresh = minY + 0.4 * (maxY - minY);
    z_inferior_thresh = minZ + 0.4 * (maxZ - minZ);

    % inferior_mask = z >= z_inferior_thresh;
    inferior_mask = z <= z_inferior_thresh;
    posterior_mask = y <= posterior_thresh;

    if strcmp(side, 'left')
        % region_mask = x <= x_left_thresh & inferior_mask & posterior_mask;
        region_mask = x <= x_left_thresh & inferior_mask & posterior_mask;
    else
        % region_mask = x >= x_right_thresh & inferior_mask & posterior_mask;
        region_mask = x >= x_right_thresh & inferior_mask & posterior_mask;
    end

    subset = branchList(region_mask, :);
    segments = unique(subset(:,4));
    best_id = NaN;
    best_length = 0;
    z_std_thresh = 3.5;

    for i = 1:length(segments)
        seg_id = segments(i);
        seg_points = subset(subset(:,4) == seg_id, :);
        num_pts = size(seg_points,1);

        if num_pts <= 40
            z_std = std(seg_points(:,3));
            if z_std < z_std_thresh && num_pts > best_length
                best_length = num_pts;
                best_id = seg_id;
            end
        else
            quarter = floor(num_pts / 4);
            parts = {seg_points(1:quarter,:); seg_points(quarter+1:2*quarter,:); ...
                     seg_points(2*quarter+1:3*quarter,:); seg_points(3*quarter+1:end,:)};

            for part = parts
                p = part{1};
                if strcmp(side, 'left')
                    x_ok = p(:,1) <= x_left_thresh;
                else
                    x_ok = p(:,1) >= x_right_thresh;
                end
                y_ok = p(:,2) <= posterior_thresh;
                % z_ok = p(:,3) >= z_inferior_thresh;
                z_ok = p(:,3) <= z_inferior_thresh;
                region_mask = x_ok & y_ok & z_ok;

                if mean(region_mask) >= 0.9
                    z_std = std(p(:,3));
                    len = size(p,1);
                    if z_std < z_std_thresh && len > best_length
                        best_length = len;
                        best_id = seg_id;
                    end
                end
            end
        end
    end
    segment_id = branchList(:,4);  % Define it here
    seg_points = branchList(segment_id == best_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);
end

function loc = extractSTRV(segment_ids, data_struct, sssv_segment_id)
    if nargin < 3
        sssv_segment_id = [];
    end
    
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);

    % Only consider points inside the 4D flow segmentation mask
    inMask = isInsideSegmentMask(data_struct, branchList);
    branchList = branchList(inMask, :);
    if isempty(branchList), loc = []; return; end

    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); segment_id = branchList(:,4);

    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);
    posterior_thresh = minY + 0.5 * (maxY - minY);
    superior_thresh  = minZ + 0.33 * (maxZ - minZ);

    posterior_mask = y <= posterior_thresh;
    % superior_mask  = z <= superior_thresh;
    superior_mask  = z >= superior_thresh;
    straight_mask = posterior_mask & superior_mask;

    straight_subset = branchList(straight_mask, :);
    straight_segments = unique(straight_subset(:,4));

    % Exclude SSSV segment if it was already assigned
    if ~isempty(sssv_segment_id) && ismember(sssv_segment_id, straight_segments)
        straight_segments = straight_segments(straight_segments ~= sssv_segment_id);
    end

    direction_threshold = 0.85;
    min_points = 15;
    % expected = [0; 1; -1]; expected = expected / norm(expected);
    expected = [0; 1; 1]; expected = expected / norm(expected);
    best_score = -Inf;
    straight_sinus_id = NaN;

    for i = 1:length(straight_segments)
        seg_id = straight_segments(i);
        if seg_id == 1, continue; end
        seg_points = straight_subset(straight_subset(:,4) == seg_id, 1:3);
        if size(seg_points,1) < min_points, continue; end

        centered = seg_points - mean(seg_points);
        [~, ~, V] = svd(centered, 'econ');
        direction = V(:,1);
        alignment = abs(dot(direction, expected));

        if alignment > direction_threshold && alignment > best_score
            best_score = alignment;
            straight_sinus_id = seg_id;
        end
    end

    if isnan(straight_sinus_id)
        loc = [];
        return;
    end

    seg_points = branchList(segment_id == straight_sinus_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);
end

function mask = isInsideSegmentMask(data_struct, points)
    % Check which centerline points lie inside the 4D flow segmentation mask.
    % points: Nx5 branchList rows [y, x, z, seg_id, cl_idx] (or Nx3 [y, x, z]).
    % Returns logical Nx1 mask: true if data_struct.segment(round(y), round(x), round(z)) > 0.
    segVol = data_struct.segment;
    dims = size(segVol);
    n = size(points, 1);
    mask = false(n, 1);
    for i = 1:n
        iy = round(points(i, 1));
        ix = round(points(i, 2));
        iz = round(points(i, 3));
        if iy >= 1 && iy <= dims(1) && ix >= 1 && ix <= dims(2) && iz >= 1 && iz <= dims(3)
            mask(i) = segVol(iy, ix, iz) > 0;
        end
    end
end

function circularity = getCircularity(data_struct, rowIdx)
    % Get circularity value for a given branchList row index
    % Returns diam_val (Rin^2/Rout^2) which is already computed
    % Value of 1 = perfect circle, <1 = irregular shape
    %
    % Inputs:
    %   data_struct - Data structure containing diam_val
    %   rowIdx - Row index in branchList (1-based)
    %
    % Output:
    %   circularity - Circularity metric (0 if invalid)
    
    if nargin < 2 || isempty(data_struct) || ~isfield(data_struct, 'diam_val')
        circularity = 0;
        return;
    end
    
    if rowIdx > 0 && rowIdx <= length(data_struct.diam_val)
        circularity = data_struct.diam_val(rowIdx);
        % Handle invalid values
        if ~isfinite(circularity) || circularity <= 0
            circularity = 0;
        end
    else
        circularity = 0;
    end
end

function bestLOC = selectICAsBASILOC(candidates, data_struct)
    % Select the best measuring point (LOC) for ICA/BASI from candidate centerline points.
    %
    % Strategy (ordered by priority):
    %   1. EXCLUDE siphon region: ICAs curve upward at their distal end (high Z). Detect the
    %      siphon onset as the point where cumulative Z starts increasing monotonically after
    %      the global Z minimum along the segment. Exclude candidates above that threshold.
    %   2. EXCLUDE area outliers: remove candidates with cross-sectional area outside
    %      [median - 2*MAD, median + 2*MAD] or below a hard floor (e.g. 0.01 cm^2).
    %   3. SCORE remaining candidates by:
    %        - Circularity (diam_val, Rin^2/Rout^2; 1 = perfect circle; higher is better)
    %        - Velocity magnitude (velMean_val; higher = better SNR / more reliable measurement)
    %      Both normalized 0–1 and combined with equal weight.
    %      Tie-break: prefer the candidate closest to the segment midpoint (avoids edges).
    %
    % Inputs:
    %   candidates  - Nx5 matrix of branchList rows [y, x, z, segment_id, point_index]
    %   data_struct - Structure with branchList, area_val, diam_val, velMean_val
    %
    % Output:
    %   bestLOC - 1x5 best candidate row, or empty if no candidates

    if isempty(candidates) || size(candidates, 2) < 5
        bestLOC = [];
        return;
    end

    branchList = data_struct.branchList;

    %% --- Map each candidate to its branchList row index ---
    nCand = size(candidates, 1);
    rowIdxs = zeros(nCand, 1);
    for i = 1:nCand
        rowIdxs(i) = find(branchList(:,4) == candidates(i,4) & ...
                          branchList(:,5) == candidates(i,5), 1);
    end
    valid = rowIdxs > 0;
    candidates = candidates(valid, :);
    rowIdxs = rowIdxs(valid);
    nCand = size(candidates, 1);
    if nCand == 0, bestLOC = []; 
        fprintf('[selectICAsBASILOC] No valid candidates found\n');
        return; 
    end

    %% --- Stage 1: Siphon exclusion ---
    % For each segment, find the Z profile and detect the siphon onset:
    % the ICA siphon is a U/S-bend where Z rises at the distal end.
    % We exclude points in the top 70th percentile of Z within their segment,
    % measured relative to the segment's own Z range.
    siphonKeep = true(nCand, 1);
    segIDs = unique(candidates(:, 4));
    allSegZ = [];
    for s = 1:numel(segIDs)
        segMask = branchList(:, 4) == segIDs(s);
        allSegZ = [allSegZ; branchList(segMask, 3)];
    end

    siphonKeep = true(nCand, 1);
    if numel(allSegZ) >= 3
        zMin = min(allSegZ);
        zMax = max(allSegZ);
        zRange = zMax - zMin;
        if zRange > eps
            siphonZThresh = zMin + 0.7 * zRange;
            siphonKeep = candidates(:, 3) <= siphonZThresh;
        end
    end

    % If siphon exclusion removes everything, relax: keep all
    if ~any(siphonKeep)
        siphonKeep = true(nCand, 1);
    end
    candidates = candidates(siphonKeep, :);
    rowIdxs = rowIdxs(siphonKeep);
    nCand = size(candidates, 1);
    if nCand == 0, bestLOC = []; fprintf('[selectICAsBASILOC] No valid candidates after siphon exclusion\n'); return; end

    %% --- Stage 2: Area outlier exclusion ---
    areaFloor = 0.01;  % cm^2 hard minimum
    areas = data_struct.area_val(rowIdxs);
    areas(~isfinite(areas)) = 0;

    medArea = median(areas(areas > 0));
    madArea = mad(areas(areas > 0), 1);  % median absolute deviation
    if isempty(medArea) || madArea == 0
        areaKeep = areas >= areaFloor;
    else
        areaKeep = areas >= max(areaFloor, medArea - 4*madArea) & ...
                   areas <= medArea + 4*madArea;
    end
    if ~any(areaKeep), areaKeep = true(nCand, 1); end  % relax if all excluded
    candidates = candidates(areaKeep, :);
    rowIdxs = rowIdxs(areaKeep);
    nCand = size(candidates, 1);
    if nCand == 0, bestLOC = []; fprintf('[selectICAsBASILOC] No valid candidates after area outlier exclusion\n'); return; end

    %% --- Stage 3: Score by circularity + flow magnitude ---
    circ = data_struct.diam_val(rowIdxs);
    % circ = computeAxialCircularity(data_struct, candidates, rowIdxs);

    circ(~isfinite(circ) | circ <= 0) = 0;
    flow = abs(data_struct.flowPerHeartCycle_val(rowIdxs));
    flow(~isfinite(flow)) = 0;
    % Normalize both to [0, 1]
    maxCirc = max(circ);
    if maxCirc > 0, normCirc = circ / maxCirc; else, normCirc = zeros(nCand,1); end
    maxFlow = max(flow);
    if maxFlow > 0, normFlow = flow / maxFlow; else, normFlow = zeros(nCand,1); end
    % Weighted combination
    wCirc = 0.4;
    wFlow = 0.6;
    score = wCirc * normCirc + wFlow * normFlow;
    [~, bestIdx] = max(score);
    bestLOC = candidates(bestIdx, :);
    fprintf('[selectICAsBASILOC] Best candidate: %f, %f, %f, %d, %d\n', bestLOC(1), bestLOC(2), bestLOC(3), bestLOC(4), bestLOC(5));
end

function circ = computeAxialCircularity(data_struct, candidates, rowIdxs)
    % Compute circularity from the Z-oriented (axial) slice of the binary mask
    % at each candidate point, instead of the reoriented cross-sectional plane.
    segVol = data_struct.segment;  % 3D binary mask [dim1 x dim2 x dim3]
    branchList = data_struct.branchList;
    nCand = numel(rowIdxs);
    circ = zeros(nCand, 1);

    for i = 1:nCand
        coord = round(branchList(rowIdxs(i), 1:3));  % [y, x, z]
        iy = coord(1); ix = coord(2); iz = coord(3);

        % Bounds check
        if iz < 1 || iz > size(segVol, 3) || ...
           iy < 1 || iy > size(segVol, 1) || ix < 1 || ix > size(segVol, 2)
            circ(i) = 0;
            continue;
        end

        % Extract axial slice at this Z
        axSlice = segVol(:, :, iz) > 0;
        if ~any(axSlice(:))
            circ(i) = 0;
            continue;
        end

        % Label connected components and pick the one closest to (y, x)
        [L, nLabels] = bwlabel(axSlice);
        if nLabels == 0
            circ(i) = 0;
            continue;
        end
        s = regionprops(L, 'Centroid');
        centroids = reshape([s.Centroid], 2, [])';  % [col, row] per component
        dists = sqrt((centroids(:,2) - iy).^2 + (centroids(:,1) - ix).^2);
        [~, closest] = min(dists);
        component = (L == closest);

        % Compute Rin^2 / Rout^2 (same as segment_cross_section_thresh)
        D = bwdist(~component);
        Rin = max(D(:));
        [xLoc, yLoc] = find(bwperim(component));
        if isempty(xLoc) || Rin == 0
            circ(i) = 0;
            continue;
        end
        Dperi = pdist2([xLoc, yLoc], [xLoc, yLoc]);
        Rout = max(Dperi(:)) / 2;
        if Rout == 0
            circ(i) = 0;
        else
            circ(i) = Rin^2 / Rout^2;
        end
    end
end

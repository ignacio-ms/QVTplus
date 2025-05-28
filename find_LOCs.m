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
        count_LICA = sum(info_LICA(:, 3) == z);
        count_RICA = sum(info_RICA(:, 3) == z);
        count_BA = sum(info_BA(:, 3) == z);

        if count_LICA == 1 && count_RICA == 1 && count_BA == 1
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
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments
    x = branchList(:,1);
    y = branchList(:,2);
    z = branchList(:,3);
    segment_id = branchList(:,4);

    dims = size(data_struct.segment);
    midY = dims(2)/3;
    midZ = dims(3)/3;

    posterior_mask = y < midY;
    superior_mask  = z < midZ;

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
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments

    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); val = branchList(:,4);

    minX = min(x); maxX = max(x);
    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);

    x_left_thresh  = minX + 0.5 * (maxX - minX);
    x_right_thresh = maxX - 0.33 * (maxX - minX);
    posterior_thresh = minY + 0.33 * (maxY - minY);
    z_inferior_thresh = minZ + 0.60 * (maxZ - minZ);

    inferior_mask = z >= z_inferior_thresh;
    posterior_mask = y <= posterior_thresh;

    if strcmp(side, 'left')
        region_mask = x <= x_left_thresh & inferior_mask & posterior_mask;
    else
        region_mask = x >= x_right_thresh & inferior_mask & posterior_mask;
    end

    subset = branchList(region_mask, :);
    segments = unique(subset(:,4));
    best_id = NaN;
    best_length = 0;
    z_std_thresh = 3.0;

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
                z_ok = p(:,3) >= z_inferior_thresh;
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

function loc = extractSTRV(segment_ids, data_struct)
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments
    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); segment_id = branchList(:,4);

    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);
    posterior_thresh = minY + 0.5 * (maxY - minY);
    superior_thresh  = minZ + 0.5 * (maxZ - minZ);

    posterior_mask = y <= posterior_thresh;
    superior_mask  = z <= superior_thresh;
    straight_mask = posterior_mask & superior_mask;

    straight_subset = branchList(straight_mask, :);
    straight_segments = unique(straight_subset(:,4));

    direction_threshold = 0.90;
    min_points = 20;
    expected = [0; 1; -1]; expected = expected / norm(expected);
    best_score = -Inf;
    straight_sinus_id = NaN;

    for i = 1:length(straight_segments)
        seg_id = straight_segments(i);
        if seg_id == 1, continue; end
        seg_points = straight_subset(straight_subset(:,4) == seg_id, 1:3);
        if size(branchList(branchList(:,4) == seg_id),1) < min_points, continue; end

        centered = seg_points - mean(seg_points);
        [~, ~, V] = svd(centered, 'econ');
        direction = V(:,1);
        alignment = abs(dot(direction, expected));

        if alignment > direction_threshold && alignment > best_score
            best_score = alignment;
            straight_sinus_id = seg_id;
        end
    end

    seg_points = branchList(segment_id == straight_sinus_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);
end

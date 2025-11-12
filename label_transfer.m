function varargout = label_transfer(funcName, varargin)
    % Dispatcher function
    switch funcName
        case 'processNeighbours'
            varargout{1} = processNeighbours(varargin{:});
        case 'assignMajorityLabels'
            varargout{1} = assignMajorityLabels(varargin{:});
        case 'expandLabels'
            varargout{1} = expandLabels(varargin{:});
        otherwise
            error('Unknown function name: %s', funcName);
    end
end


function validLabels = processNeighbours(nearestIdx, dist, multiLabels, threshold, k)
    % Find valid neighbors within distance threshold
    validIdx = nearestIdx(dist <= threshold);
    validLabels = multiLabels(validIdx);
    
    % Pad with NaNs if fewer than k neighbors
    if length(validLabels) < k
        validLabels = [validLabels; NaN(k - length(validLabels), 1)];
    else
        validLabels = validLabels(1:k); % Trim to k neighbors
    end
end

function majorityLabels = assignMajorityLabels(validNeighbors, segMatrix, rows, cols, slices)
    % Assign the majority label based on valid neighbors
    majorityLabels = zeros(length(validNeighbors), 1);
    for i = 1:length(validNeighbors)
        neighborLabels = validNeighbors{i}(~isnan(validNeighbors{i}) & validNeighbors{i} ~= 0);
        if ~isempty(neighborLabels)
            majorityLabels(i) = mode(neighborLabels);
        else
            % If no valid multilabel segments, assign 100 * current label
            majorityLabels(i) = segMatrix(rows(i), cols(i), slices(i)) * 100;
        end
    end
end

function updatedMatrix = expandLabels(matrix, targetLabel, threshold)
    % Expand the labels for a given target label
    [rows, cols, slices] = ind2sub(size(matrix), find(matrix == targetLabel));
    points = [rows, cols, slices];
    
    [nonzeroRows, nonzeroCols, nonzeroSlices] = ind2sub(size(matrix), find(matrix > 0 & matrix ~= targetLabel));
    nonzeroPoints = [nonzeroRows, nonzeroCols, nonzeroSlices];
    nonzeroLabels = matrix(matrix > 0 & matrix ~= targetLabel);
    
    % Perform k-NN search
    [nearestIndices, distances] = knnsearch(nonzeroPoints, points, 'K', 10);
    
    % Assign new labels
    newLabels = zeros(size(points, 1), 1);
    for i = 1:size(points, 1)
        validIndices = nearestIndices(i, distances(i, :) <= threshold);
        validLabels = nonzeroLabels(validIndices);
        if ~isempty(validLabels)
            newLabels(i) = mode(validLabels);
        else
            newLabels(i) = targetLabel; % Retain original label if no valid neighbors
        end
    end
    
    % Update the matrix
    linearIndices = sub2ind(size(matrix), rows, cols, slices);
    updatedMatrix = matrix;
    updatedMatrix(linearIndices) = newLabels;
end

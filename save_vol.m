function save_vol(output_path, data_struct, LOCs)
    % Load reference image
    imageData = spm_vol(fullfile(output_path, 'QVT_MAG.nii'));
    fullSize = imageData.dim;
    volume = zeros(fullSize);

    % Initialize branchList with -1 labels
    branchList = data_struct.branchList;
    branchList(:,4) = -1;

    % Define label mapping
    label_map = struct( ...
        'BASI',   1, ...
        'LICA', 2, ...
        'RICA', 3, ...
        'LMCA', 4, ...
        'RMCA', 5, ...
        'LACA', 6, ...
        'RACA', 7, ...
        'LPCA', 8, ...
        'RPCA', 9, ...
        'SSSV', 10, ...
        'STRV', 11, ...
        'LTSV', 12, ...
        'RTSV', 13 ...
    );

    % Step 1: Collect branch usage info
    fields = fieldnames(LOCs);
    branch_to_vessels = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for i = 1:length(fields)
        field = fields{i};
        info = LOCs.(field);
        branch_id = info(1);
        seed_idx = info(2);
    
        if ~isKey(branch_to_vessels, branch_id)
            branch_to_vessels(branch_id) = {};
        end
    
        vessel_list = branch_to_vessels(branch_id);
        vessel_list{end+1} = struct( ...
            'field', field, ...
            'label', label_map.(field), ...
            'seed_idx', seed_idx ...
        );
        branch_to_vessels(branch_id) = vessel_list;
    end

    % Step 2: Label each branch
    original_ids = data_struct.branchList(:,4);
    coords = data_struct.branchList(:,1:3);

    for branch_id = branch_to_vessels.keys
        vessel_infos = branch_to_vessels(branch_id{1});
        idx_in_branch = find(round(original_ids) == branch_id{1});

        if length(vessel_infos) == 1
            % Unique: assign entire branch to that label
            label = vessel_infos{1}.label;
            branchList(idx_in_branch, 4) = label;

        elseif length(vessel_infos) == 2
            % Shared branch: split via distance from seeds
            v1 = vessel_infos{1};
            v2 = vessel_infos{2};

            % Get coordinates of the full branch
            branch_coords = coords(idx_in_branch, :);
            seed1_coord = coords(v1.seed_idx, :);
            seed2_coord = coords(v2.seed_idx, :);

            % Compute distances
            d1 = sum((branch_coords - seed1_coord).^2, 2);
            d2 = sum((branch_coords - seed2_coord).^2, 2);

            % Assign based on closest seed
            for i = 1:length(idx_in_branch)
                row = idx_in_branch(i);
                if d1(i) <= d2(i)
                    branchList(row, 4) = v1.label;
                else
                    branchList(row, 4) = v2.label;
                end
            end
        else
            warning('Branch %d is shared by more than 2 vessels — not supported.', branch_id{1});
        end
    end

    % Step 3: Write to volume
    x = round(branchList(:,1));
    y = round(branchList(:,2));
    z = round(branchList(:,3));
    val = branchList(:,4);

    for i = 1:length(x)
        xi = x(i); yi = y(i); zi = z(i);
        if xi > 0 && yi > 0 && zi > 0 && ...
           xi <= fullSize(1) && yi <= fullSize(2) && zi <= fullSize(3)
            volume(xi, yi, zi) = val(i);
        end
    end

    % Save as NIfTI
    V = imageData;
    V.fname = fullfile(output_path, 'branchmask.nii');
    V.descrip = 'Multilabel mask with split shared branches';
    spm_write_vol(V, volume);
end

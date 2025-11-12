%save QVT image and segmentation
MAG = imageData.MAG;
seg = imageData.Segmented;

% Write the volume to a NIfTI file


% Assume V is the volume you want to save, and fname is the filename
V = struct;
V.fname = 'QVT_mag.nii';  % Set the filename in the NIfTI header
V.dim = size(MAG);  % Set the dimensions of the volume
V.dt = [spm_type('float32'), 0];  % Set the data type of the volume
V.mat = eye(4);  % Set the affine transformation matrix
V.descrip = 'NIfTI volume';  % Set a description
spm_write_vol(V, MAG);  % Write the volume to a NIfTI file

V = struct;
V.fname = 'QVT_seg.nii';  % Set the filename in the NIfTI header
V.dim = size(seg);  % Set the dimensions of the volume
V.dt = [spm_type('uint8'), 0];  % Set the data type of the volume
V.mat = eye(4);  % Set the affine transformation matrix
V.descrip = 'NIfTI volume';  % Set a description
spm_write_vol(V, seg);  % Write the volume to a NIfTI file


% Load the NIfTI file
NII = spm_vol('BMRI198894_4DQflowNeuro_AP_1803_tAvgProj_resampled.nii');
NII = spm_vol('BMRI198894_4DQflowNeuro_AP_1803_tAvgProj_eICAB_WB.nii');

NII = spm_vol('eICAB_seg_realign.nii');


% Read the volume
volume = spm_read_vols(NII);
img2= flip(volume,3);

dims=size(img2);

A_rotated = zeros(dims(2),dims(1),dims(3));  % Preallocate the rotated matrix

% Loop over each slice along the z-dimension
for k = 1:size(img2, 3)
    % Rotate the k-th slice 90 degrees clockwise
    A_rotated(:, :, k) = rot90(img2(:, :, k), 1);
end

seg2 = flip(A_rotated,2);
NII.dim = dims([2 1 3]);
new_voxel_size = [1 1 1];  % Replace with your new voxel size
NII.mat(1:3, 1:3) = diag(new_voxel_size);

% Change the origin
new_origin = [0 0 0];  % Replace with your new origin
NII.mat(1:3, 4) = new_origin;
NII.fname = 'eICAB_seg.nii';

spm_write_vol(NII, seg2);

% realign 4DFlow
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.realign.estwrite.data = {{'eICAB_vol.nii'},{'eICAB_vol'}};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {'QVT_mag.nii'};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch


% coregister
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'QVT_mag.nii'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'rBMRI198894_4DQflowNeuro_AP_1803_tAvgProj_resampled.nii'};%{'eICAB_vol.nii'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {'rBMRI198894_4DQflowNeuro_AP_1803_tAvgProj_eICAB_WB.nii'};%{'eICAB_seg.nii'};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);


newseg2 = spm_vol('reICAB_seg.nii');
qvtseg = spm_vol('QVT_seg.nii');

finalseg = spm_vol('rrBMRI198894_4DQflowNeuro_AP_1803_tAvgProj_eICAB_WB.nii');

volshow(double(finalseg.private.dat).*double(qvtseg.private.dat))
%i can read and save the vols through spm

orig_seg = imageData.Segmented;
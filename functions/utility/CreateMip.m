% This function just makes a maximum intensity projection of Time of Flight
% (or any other white blood MRA). Put in directory of the folder. 
%TODO: Create 4D Flow MIP on running ParamMap.
Path2TOF='';
DIR=dir(Path2Mip);
filename=DIR(3).name;
slice = dicomread(fullfile(Path2Mip,filename));
MIP=zeros([size(slice) 2]);
for i=3:length(DIR)
    filename=DIR(i).name;
    slice = dicomread(fullfile(Path2Mip,filename));
    MIP(:,:,2)=slice;
    MIP(:,:,1)=max(MIP,[],3);
end
imagesc(MIP(:,:,1))
colormap(gray)
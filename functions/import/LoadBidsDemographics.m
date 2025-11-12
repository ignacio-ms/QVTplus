% This function will load demographics for your study cohort if they're
% stored in a BIDS formatted sourcedata directory. It reads the first DICOM
% file for each subject and collects demographics from DICOM header. 

%Note: The organisation is based on alphabetical file order. If you want to
%reorganise another way, add functionality.
%% Initialization
path2bids='C:\Users\sdem348\Desktop\Dempsey2023MultiRes_Cohort';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Don't change below %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Or do
path2subs=fullfile(path2bids,'sourcedata');
DIR=dir(path2subs);
Demos=cell([(length(DIR)-2),4]);
for i=3:length(DIR(:))
    path2patient=fullfile(path2subs,DIR(i).name);
    DIR2=dir(path2patient);
    DIR3=dir(fullfile(path2patient,DIR2(3).name));
    filename=DIR3(3).name;
    INFO=dicominfo(fullfile(fullfile(path2patient,DIR2(3).name),filename));
    Demos(i-2,1)={DIR(i).name};
    Demos(i-2,2)={INFO.PatientAge(2:3)};
    Demos(i-2,3)={INFO.PatientSex};
    Demos(i-2,4)={INFO.PatientWeight};
    Demos(i-2,5)={INFO.HeartRate};
    Demos(i-2,6)={INFO.HeartRate};
    Demos(i-2,7)={INFO.HeartRate};
    Demos(i-2,8)={INFO.HeartRate};
    Demos(i-2,9)={INFO.HeartRate};
    Demos(i-2,10)={INFO.HeartRate};
end
save(fullfile(path2bids,'derivatives','QVT','population','Demographics.mat'),'Demos')
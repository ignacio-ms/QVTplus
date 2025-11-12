clc;clear;
%This function will rename exported DICOM folders by their scan name in the
%DICOM metadata. Good for visualising during folder exploration, but
%unecessary to any processing. 
path2data='C:\Users\sdem348\Desktop\Shaihu';
DIR=dir(path2data);
EXP='s\d*';
for i=8:length(DIR)
    foldername1=DIR(i).name;
    if regexp(foldername1,EXP) == 1
        DIR2=dir(fullfile(path2data,foldername1));
        path2dcms=fullfile(path2data,foldername1,DIR2(3).name,foldername1);
        DIR3=dir(path2dcms);
        try
            INFO=dicominfo(fullfile(path2dcms,DIR3(3).name));
            Name=INFO.SeriesDescription;
            Name = erase(Name,":");
                try
                    movefile(path2dcms,fullfile(path2data,Name))
                    catch
                end
        catch
        end
    end
end

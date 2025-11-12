function [] = autoCollectFlow(SVPATH)
% autoCollectFlow loads an interactive table to initialise any locations or
% vessels to process for PITC or PI pipelines. It replaces the picking and
% saving of standard QVT (which takes time to save points into excel and
% then demands being read from that excel). 
%
% This way, we have much more control on the processing, and can pull right
% from the saved data_struct
%
% Used by: paramMap.m
% Dependencies: readLabels.m
    Artery = ["L ICA";"R ICA";"L MCA";"R MCA";"L ACA";"R ACA";"L PCA";"R PCA";"BA";"Exclude"];
    Loc = [" ";" ";" ";" ";" ";" ";" ";" ";" ";" "];
    P = table('Size',[10,3],'VariableTypes',{'cellstr','cellstr','cellstr'},'VariableNames',["Artery","Label","Loc"]);
    P.Artery=Artery;
    P.Label=Loc;
    P.Loc=Loc;
    if exist(fullfile(SVPATH,'LabelsQVT.csv')) == 2
        Values = readLabels(SVPATH);
        P.Label=Values(:,2);
        P.Loc=Values(:,3);
    end
    % fig = uifigure("Position",[100 100 350 375]);
    % fig.Tag = 'testGUI2_tag';
    % uit = uitable(fig,"Data",P);
    % uit.ColumnEditable=true;
    % But = uibutton(fig,"Text","Done","Position",[150 340 50 25],...
    %     "ButtonPushedFcn", @(src,event) plotButtonPushed(fig,uit));
    % function plotButtonPushed(fig,uit)
    %     writetable(uit.Data,fullfile(SVPATH,'LabelsQVT.csv'));
    %     close(fig)
    % end
end
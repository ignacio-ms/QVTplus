function [PWV,R]=enc_PWV(data_struct,PI_scat,time,Labels,params)
    PWV=zeros([1 12]);
    R=zeros([1 12]);
    Flows=data_struct.flowPulsatile_val;
    if params.PWV(1)==1
        %=======================================
        %Process Bjornfoot 2021 Individual Roots Bjornfoot Weights
        %=======================================
        for root=1:3
            Vals=squeeze(PI_scat(:,[1 3 4 5],root));
            [row,~]=find(Vals(:,2)==0);
            Vals(row,:)=[];
            if ~isempty(Vals(:,1))
                [pwv,A]=enc_PWVBjornfoot(Flows,Vals,time,0,params);
                PWV(1,root) = pwv;
                W{root}=A;
            else
                PWV(1,root) = pwv;
            end
        end
    end
    if params.PWV(2)==1
        %=====================================
        %Process Fielding (Xcor) Root specific  Bjornfoot Weights
        %=====================================
        if params.PltFlag==1
            figure()
        end
        for root=1:3
            Vals=squeeze(PI_scat(:,[1 3 4 5],root)); %D, ID, Q, A
            [row,~]=find(Vals(:,2)==0);
            Vals(row,:)=[];
            if ~isempty(Vals(:,1))
                [pwv,r]=enc_PWV_XCor(Vals,data_struct,Labels,time,0,W{root},params,root);
                PWV(1,root+3) = pwv;
                R(1,root+3) = r;
            else
                PWV(1,root+3) = 0;
            end
        end
        if params.PltFlag==1
            if params.SaveData==1
                path2data=string(fullfile(params.data_dir,params.subject));
                saveas(gcf,fullfile(path2data,'PWV_bw.jpg'))
                close(gcf)
            end    
        end
    end

    if params.PWV(3)==1
        %=======================================
        %Process Bjornfoot 2021 Individual Roots, Dempsey Weights
        %=======================================
        for root=1:3
            Vals=squeeze(PI_scat(:,[1 3 4 5],root));
            [row,~]=find(Vals(:,2)==0);
            Vals(row,:)=[];
            if ~isempty(Vals(:,1))
                [pwv,~]=enc_PWVBjornfoot(Flows,Vals,time,1,params);
                PWV(1,root+6) = pwv;
            else
                PWV(1,root+6) = pwv;
            end
        end
    end

    if params.PWV(4)==1
        %================================================
        %Process Fielding (Xcor) Root specific (all Data) Dempsey Weights
        %================================================
        if params.PltFlag==1
            figure()
        end
        for root=1:3
            Vals=squeeze(PI_scat(:,[1 3 4 5],root)); %D, ID, Q, A
            [row,~]=find(Vals(:,2)==0);
            Vals(row,:)=[];
            if ~isempty(Vals(:,1))
                [pwv,r]=enc_PWV_XCor(Vals,data_struct,Labels,time,1,[],params,root);
                PWV(1,root+9) = pwv;
                R(1,root+9) = r;
            else
                PWV(1,root+9) = 0;
            end
        end
        if params.PltFlag==1
            if params.SaveData==1
                path2data=string(fullfile(params.data_dir,params.subject));
                saveas(gcf,fullfile(path2data,'PWV_dw.jpg'))
                close(gcf)
            end    
        end          
    end

end

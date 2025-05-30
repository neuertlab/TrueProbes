function [MapF_DPS, MapF_SiteLoc,MapF_Vars,Num_of_Molecule_Sites] = combineMultipleBindingSiteMaps_V2(Map_DPS, Map_SiteLoc, Map_Vars)
% Check input validity
%        MapMultiVars{k}{1} = updated_Multi_KbMod{k};
%        MapMultiVars{k}{2} = updated_Multi_TmMod{k};
%        MapMultiVars{k}{3} = updated_Multi_dHeqMod{k};
%        MapMultiVars{k}{4} = updated_Multi_dSeqMod{k};
%        MapMultiVars{k}{5} = updated_Multi_dCpMod{k};
%        MapMultiVars{k}{6} = updated_Multi_dHfMod{k};
%        MapMultiVars{k}{7} = updated_Multi_dSfMod{k};
%        MapMultiVars{k}{8} = updated_Multi_dHrMod{k};
%        MapMultiVars{k}{9} = updated_Multi_dSrMod{k}; 
%        MapMultiVars{k}{10} = updated_Multi_KbComp{k};
%        MapMultiVars{k}{11} = updated_Multi_dHeq_Complement{k};
%        MapMultiVars{k}{12} = updated_Multi_dSeq_Complement{k};
%        MapMultiVars{k}{13} = updated_Multi_dCp_Complement{k};
%        MapMultiVars{k}{14} = updated_Multi_dHf_Complement{k};
%        MapMultiVars{k}{15} = updated_Multi_dSf_Complement{k};
%        MapMultiVars{k}{16} = updated_Multi_dHr_Complement{k};
%        MapMultiVars{k}{17} = updated_Multi_dSr_Complement{k};
if ~iscell(Map_DPS) || ~iscell(Map_SiteLoc) || length(Map_DPS) ~= length(Map_SiteLoc)
    error('Input should be cell arrays of equal length');
end
%Correlary function that given new and old siteLoc can tell you Sold<->Snew
% Get number of maps
num_maps = length(Map_DPS);
max_map_sites = cellfun(@(x) size(x,3),Map_DPS);
Num_of_Molecule_Sites = zeros(1,size(Map_DPS{1},2));

[~,ord] = sort(max_map_sites,'descend');
Map_DPS = Map_DPS(ord);
Map_SiteLoc = Map_SiteLoc(ord);
Map_Vars = Map_Vars(ord);

% Initialize Map3_DPS and Map3_SiteLoc
MapQ_DPS = Map_DPS{1};
MapQ_SiteLoc = Map_SiteLoc{1};
MapQ_Vars = Map_Vars{1};
%MapQ_Vars{k} =  ndSparse.build([size(MapQ_SiteLoc,1),size(MapQ_SiteLoc,2),size(MapQ_SiteLoc,3)],0);%P T S
MapCurr_DPS = Map_DPS{1};
MapCurr_SiteLoc = Map_SiteLoc{1};
MapCurr_Vars = Map_Vars{1};       
%by difinition site cutting is with overlap regions, to get max and min
%for probes involved in those sites it is accurate since siteLoc references
%each map and probes individually. which can resolve any new redundancies

    %NEW SITE MAP
    %Probes are only in Map 1 bind some site X
    %Probes are only in Map 2 bind some site X
    %Probes are in Map1 & 2 bind some site X
    %[PTSi]->[PTSf]
        %[P,S] -> [P',S'] linIndex ->linIndex' fix cases for KbComp when
        %map_idx is first value or not for assigning new values
for map_idx = 2:num_maps
    MapProp_DPS = Map_DPS{map_idx};
    MapProp_SiteLoc = Map_SiteLoc{map_idx};% Combine the site locations from Map1 and Map2
    MapProp_Vars = Map_Vars{map_idx};
    PropLinIdx = find(MapProp_DPS);
    CurrLinIdx = find(MapCurr_DPS);
    %[E Pi Ti Si Ax Bx Li Az Bz Zi] Z new sites
    [probeProp,targetProp,SiteProp] = ind2sub(size(MapProp_DPS),PropLinIdx);
    Em_Prop = map_idx*ones(length(probeProp),1);
    Ax_Prop = full(MapProp_SiteLoc(sub2ind(size(MapProp_SiteLoc),probeProp,targetProp,SiteProp,ones(1,length(targetProp))')));
    Bx_Prop = full(MapProp_SiteLoc(sub2ind(size(MapProp_SiteLoc),probeProp,targetProp,SiteProp,2*ones(1,length(targetProp))')));
    Prop_Matrix = [Em_Prop probeProp targetProp SiteProp Ax_Prop Bx_Prop 0*ones(length(targetProp),3) PropLinIdx 0*ones(length(targetProp),1)];
    [probeCurr,targetCurr,SiteCurr] = ind2sub(size(MapCurr_DPS),CurrLinIdx);
    Em_Curr = (map_idx-1)*ones(length(probeCurr),1); 
    Ax_Curr = full(MapCurr_SiteLoc(sub2ind(size(MapCurr_SiteLoc),probeCurr,targetCurr,SiteCurr,ones(1,length(targetCurr))')));
    Bx_Curr = full(MapCurr_SiteLoc(sub2ind(size(MapCurr_SiteLoc),probeCurr,targetCurr,SiteCurr,2*ones(1,length(targetCurr))')));
    Curr_Matrix = [Em_Curr probeCurr targetCurr SiteCurr Ax_Curr Bx_Curr 0*ones(length(targetCurr),3) CurrLinIdx 0*ones(length(targetCurr),1)];
    targetsOnlyInProp = setdiff(targetProp,targetCurr);
    targetsOnlyInCurr = setdiff(targetCurr,targetProp);
    targetsInBoth = intersect(targetProp,targetCurr);
    
    Join_Matrix = [Prop_Matrix;Curr_Matrix];
    [~,ic,~] = unique(Join_Matrix(:,2:6),'rows','stable');
    Join_Matrix = Join_Matrix(ic,:);  
    for i = 1:length(targetsInBoth)        
        En = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),1);
        Pn = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),2);
        Tn = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),3);
        Sn = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),4);       
        An = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),5);
        Bn = Join_Matrix(Join_Matrix(:,3)==targetsInBoth(i),6);
        Cn = [unique([An Bn]).'];
        probesInInterval = cell(1,length(Cn));
        for x = 1:length(Cn)
            probesInInterval{x} = find(ge(Cn(x)-An,0).*ge(Bn-An,Cn(x)-An)).';
        end
        %Gets sites with unique sets of probes
        charIntervalArray = cellfun(@num2str, probesInInterval,'Un',0);
        [charIntervalArray_Unique,~,~] = unique(charIntervalArray,'stable');
        probesInInterval_Unique = cellfun(@str2num,charIntervalArray_Unique,'Un',0);
        %Lists number of probes binding in a site
        InInterval_order = cellfun(@length,probesInInterval_Unique);
        %Gets boundary of probe site 
        InInterval_Cmax = cell2mat(cellfun(@(x) max([An(x).' Bn(x).']),probesInInterval_Unique,'UniformOutput',false));
        InInterval_Cmin = cell2mat(cellfun(@(x) min([An(x).' Bn(x).']),probesInInterval_Unique,'UniformOutput',false));    
        IsSubSet = ndSparse.build([length(probesInInterval_Unique),length(probesInInterval_Unique)],0);
        %get proposed sites which are a subset of another
        points_belowMax = find(InInterval_order<max(InInterval_order)); %Finds points below maximum size
        for v = points_belowMax
            %subset to a larger interval
            points_req1 = find(InInterval_order>InInterval_order(v));   
            %and contained within that interval 
            %i.e. bounds of subset I within larger set J
            points_req = points_req1(find(ge(InInterval_Cmax(points_req1),InInterval_Cmax(v)).*ge(InInterval_Cmin(v),InInterval_Cmin(points_req1))));
            if (~isempty(points_req))
                IsSubSet(v,points_req) = cell2mat(cellfun(@(x) all(ismember(probesInInterval_Unique{v},x)),{probesInInterval_Unique{points_req}},'UniformOutput',false));
            end
        end            
        Sx = {probesInInterval_Unique{sum(IsSubSet,2)==0}};
        Num_of_Molecule_Sites(targetsInBoth(i)) = length(Sx);
        E_List = cell2mat(cellfun(@(x) En(x)',Sx,'Un',0))';
        P_List = cell2mat(cellfun(@(x) Pn(x)',Sx,'Un',0))';
        Ti_List = cell2mat(cellfun(@(x) Tn(x)',Sx,'Un',0))';
        Si_List = cell2mat(cellfun(@(x) Sn(x)',Sx,'Un',0))';
        Sz_List = cell2mat(arrayfun(@(x) x*ones(1,length(Sx{x})),1:length(Sx),'Un',0))';
        Ai_List = cell2mat(cellfun(@(x) An(x)',Sx,'Un',0))';
        Bi_List = cell2mat(cellfun(@(x) Bn(x)',Sx,'Un',0))';
        MatchLoc_Func = @(x) find((E_List(x)-Join_Matrix(:,1)).^2+(P_List(x)-Join_Matrix(:,2)).^2+(Ti_List(x)-Join_Matrix(:,3)).^2+(Si_List(x)-Join_Matrix(:,4)).^2 == min((E_List(x)-Join_Matrix(:,1)).^2+(P_List(x)-Join_Matrix(:,2)).^2+(Ti_List(x)-Join_Matrix(:,3)).^2+(Si_List(x)-Join_Matrix(:,4)).^2));   
        JoinMatrixLoc = arrayfun(@(x) MatchLoc_Func(x),1:length(Sz_List)); 
        Join_Matrix(JoinMatrixLoc,7) = Sz_List;
        Join_Matrix(JoinMatrixLoc,8) = Ai_List;
        Join_Matrix(JoinMatrixLoc,9) = Bi_List;
        if (max(Join_Matrix(JoinMatrixLoc,7))>size(MapQ_DPS,3))
           MapQ_DPS(1,1,max(Join_Matrix(JoinMatrixLoc,7)))=0;  
           MapQ_SiteLoc(1,targetsInBoth(i),max(Join_Matrix(JoinMatrixLoc,7)),1)=0; 
           for k = 1:9
           MapQ_Vars{k}(1,targetsInBoth(i),max(Join_Matrix(JoinMatrixLoc,7)),1)=0;  
           end     
        end
        UpdatesInProp = find(Join_Matrix(JoinMatrixLoc,1)==map_idx);
        pI_VecP = Join_Matrix(JoinMatrixLoc(UpdatesInProp),2);
        sI_VecP = Join_Matrix(JoinMatrixLoc(UpdatesInProp),4);
        sF_VecP = Join_Matrix(JoinMatrixLoc(UpdatesInProp),7);
        linDPS_Pi = sub2ind(size(MapProp_DPS),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP);
        linDPS_Po = sub2ind(size(MapQ_DPS),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP);
        MapQ_DPS(linDPS_Po) = MapProp_DPS(linDPS_Pi);     
        linSiteLocA_Po = sub2ind(size(MapQ_SiteLoc),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,1*ones(length(pI_VecP),1));
        linSiteLocB_Po = sub2ind(size(MapQ_SiteLoc),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,2*ones(length(pI_VecP),1));
        MapQ_SiteLoc(linSiteLocA_Po) = Join_Matrix(JoinMatrixLoc(UpdatesInProp),8);
        MapQ_SiteLoc(linSiteLocB_Po) = Join_Matrix(JoinMatrixLoc(UpdatesInProp),9);   
        UpdatesInCurr = find(Join_Matrix(JoinMatrixLoc,1)==map_idx-1);
        pI_VecC = Join_Matrix(JoinMatrixLoc(UpdatesInCurr),2);
        sI_VecC = Join_Matrix(JoinMatrixLoc(UpdatesInCurr),4);
        sF_VecC = Join_Matrix(JoinMatrixLoc(UpdatesInCurr),7);
        linDPS_Ci = sub2ind(size(MapCurr_DPS),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC);
        linDPS_Co = sub2ind(size(MapQ_DPS),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC);
        MapQ_DPS(linDPS_Co) = MapCurr_DPS(linDPS_Ci);   
        linSiteLocA_Co = sub2ind(size(MapQ_SiteLoc),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,1*ones(length(pI_VecC),1));
        linSiteLocB_Co = sub2ind(size(MapQ_SiteLoc),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,2*ones(length(pI_VecC),1));
        MapQ_SiteLoc(linSiteLocA_Co) = Join_Matrix(JoinMatrixLoc(UpdatesInCurr),8);
        MapQ_SiteLoc(linSiteLocB_Co) = Join_Matrix(JoinMatrixLoc(UpdatesInCurr),9);   
        for m = 1:size(MapProp_Vars{1},4)
        %KbMod 
        linVar1_Pi = sub2ind(size(MapProp_Vars{1}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar1_Po = sub2ind(size(MapQ_Vars{1}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{1}(linVar1_Po) = MapProp_Vars{1}(linVar1_Pi);
        linVar1_Ci = sub2ind(size(MapCurr_Vars{1}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar1_Co = sub2ind(size(MapQ_Vars{1}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{1}(linVar1_Co) = MapCurr_Vars{1}(linVar1_Ci);
        %TmMod
        linVar2_Pi = sub2ind(size(MapProp_Vars{2}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar2_Po = sub2ind(size(MapQ_Vars{2}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{2}(linVar2_Po) = MapProp_Vars{2}(linVar2_Pi);
        linVar2_Ci = sub2ind(size(MapCurr_Vars{2}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar2_Co = sub2ind(size(MapQ_Vars{2}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{2}(linVar2_Co) = MapCurr_Vars{2}(linVar2_Ci);     
        %dHeqMod,dSeqMod,dCpMod
        linVar3_Pi = sub2ind(size(MapProp_Vars{3}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar3_Po = sub2ind(size(MapQ_Vars{3}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{3}(linVar3_Po) = MapProp_Vars{3}(linVar3_Pi);
        linVar3_Ci = sub2ind(size(MapCurr_Vars{3}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar3_Co = sub2ind(size(MapQ_Vars{3}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{3}(linVar3_Co) = MapCurr_Vars{3}(linVar3_Ci);   
        linVar4_Pi = sub2ind(size(MapProp_Vars{4}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar4_Po = sub2ind(size(MapQ_Vars{4}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{4}(linVar4_Po) = MapProp_Vars{4}(linVar4_Pi);
        linVar4_Ci = sub2ind(size(MapCurr_Vars{4}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar4_Co = sub2ind(size(MapQ_Vars{4}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{4}(linVar4_Co) = MapCurr_Vars{4}(linVar4_Ci);   
        linVar5_Pi = sub2ind(size(MapProp_Vars{5}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar5_Po = sub2ind(size(MapQ_Vars{5}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{5}(linVar5_Po) = MapProp_Vars{5}(linVar5_Pi);
        linVar5_Ci = sub2ind(size(MapCurr_Vars{5}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar5_Co = sub2ind(size(MapQ_Vars{5}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{5}(linVar5_Co) = MapCurr_Vars{5}(linVar5_Ci);   
        %dHfMod,dSfMod,dHrMod,dSrMod
        linVar6_Pi = sub2ind(size(MapProp_Vars{6}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar6_Po = sub2ind(size(MapQ_Vars{6}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{6}(linVar6_Po) = MapProp_Vars{6}(linVar6_Pi);
        linVar6_Ci = sub2ind(size(MapCurr_Vars{6}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar6_Co = sub2ind(size(MapQ_Vars{6}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{6}(linVar6_Co) = MapCurr_Vars{6}(linVar6_Ci);  
        linVar7_Pi = sub2ind(size(MapProp_Vars{7}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar7_Po = sub2ind(size(MapQ_Vars{7}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{7}(linVar7_Po) = MapProp_Vars{7}(linVar7_Pi);
        linVar7_Ci = sub2ind(size(MapCurr_Vars{7}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar7_Co = sub2ind(size(MapQ_Vars{7}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{7}(linVar7_Co) = MapCurr_Vars{7}(linVar7_Ci);
        linVar8_Pi = sub2ind(size(MapProp_Vars{8}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar8_Po = sub2ind(size(MapQ_Vars{8}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{8}(linVar8_Po) = MapProp_Vars{8}(linVar8_Pi);
        linVar8_Ci = sub2ind(size(MapCurr_Vars{8}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar8_Co = sub2ind(size(MapQ_Vars{8}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{8}(linVar8_Co) = MapCurr_Vars{8}(linVar8_Ci);
        linVar9_Pi = sub2ind(size(MapProp_Vars{9}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP,m);
        linVar9_Po = sub2ind(size(MapQ_Vars{9}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP,m);
        MapQ_Vars{9}(linVar9_Po) = MapProp_Vars{9}(linVar9_Pi);
        linVar9_Ci = sub2ind(size(MapCurr_Vars{9}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC,m);
        linVar9_Co = sub2ind(size(MapQ_Vars{9}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC,m);
        MapQ_Vars{9}(linVar9_Co) = MapCurr_Vars{9}(linVar9_Ci);
        end    
%         linVar2_Pi = sub2ind(size(MapProp_Vars{2}),targetsInBoth(i)*ones(length(pI_VecP),1),sI_VecP);
%         linVar2_Po = sub2ind(size(MapQ_Vars{2}),pI_VecP,targetsInBoth(i)*ones(length(pI_VecP),1),sF_VecP);
%         MapQ_Vars{2}(linVar2_Po) = MapProp_Vars{2}(linVar2_Pi);  
%         if (map_idx==2)
%         linVar2_Ci = sub2ind(size(MapCurr_Vars{2}),targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC);
%         else
%         linVar2_Ci = sub2ind(size(MapCurr_Vars{2}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sI_VecC);
%         end
%         linVar2_Co = sub2ind(size(MapQ_Vars{2}),pI_VecC,targetsInBoth(i)*ones(length(pI_VecC),1),sF_VecC); 
%         MapQ_Vars{2}(linVar2_Co) = MapCurr_Vars{2}(linVar2_Ci);
        %filtering site overlap redundancy 
        for p=1:size(MapQ_DPS,1)
           Iz = find(diff([0 full(reshape(MapQ_DPS(p,targetsInBoth(i),:),[1 size(MapQ_DPS,3)])) 0])>0);
           if (~isempty(Iz))
               MapQ_DPS(p,targetsInBoth(i),setdiff(1:size(MapQ_DPS,3),Iz)) = 0;
               MapQ_Vars{1}(p,targetsInBoth(i),setdiff(1:size(MapQ_DPS,3),Iz)) = 0;      
           end
        end
    end    
    %MapQ_DPS(:,targetsOnlyInProp,:) = MapProp_DPS(:,targetsOnlyInProp,:);
    linP = find(MapProp_DPS(:,targetsOnlyInProp,:));
    [p1,t1,s1] = ind2sub(size(MapProp_DPS(:,targetsOnlyInProp,:)),linP);
    linPi = sub2ind(size(MapProp_DPS),p1,targetsOnlyInProp(t1),s1);
    linPo = sub2ind(size(MapQ_DPS),p1,targetsOnlyInProp(t1),s1);
    MapQ_DPS(linPo) = MapProp_DPS(linPi);%slow?
    %MapQ_DPS(:,targetsOnlyInCurr,1:numSitesCurr) = MapCurr_DPS(:,targetsOnlyInCurr,1:numSitesCurr); 
    linC = find(MapCurr_DPS(:,targetsOnlyInCurr,:));
    [p2,t2,s2] = ind2sub(size(MapCurr_DPS(:,targetsOnlyInCurr,:)),linC);
    linCi = sub2ind(size(MapCurr_DPS),p2,targetsOnlyInCurr(t2),s2);
    linCo = sub2ind(size(MapQ_DPS),p2,targetsOnlyInCurr(t2),s2);
    MapQ_DPS(linCo) = MapCurr_DPS(linCi);%slow?
    %MapQ_SiteLoc(:,targetsOnlyInProp,1:numSitesProp,:) = MapProp_SiteLoc(:,targetsOnlyInProp,1:numSitesProp,:);
    MapProp_SiteLocA = squeeze(MapProp_SiteLoc(:,:,:,1));
    linPA = find(MapProp_SiteLocA(:,targetsOnlyInProp,:));
    [p1,t1,s1] = ind2sub(size(MapProp_SiteLocA(:,targetsOnlyInProp,:)),linPA);
    linPAi = sub2ind(size(MapProp_SiteLoc),p1,targetsOnlyInProp(t1),s1,1*ones(length(t1),1));
    linPAo = sub2ind(size(MapQ_SiteLoc),p1,targetsOnlyInProp(t1),s1,1*ones(length(t1),1));
    MapQ_SiteLoc(linPAo) = MapProp_SiteLoc(linPAi);
    MapProp_SiteLocB = squeeze(MapProp_SiteLoc(:,:,:,2));
    linPB = find(MapProp_SiteLocB(:,targetsOnlyInProp,:));
    [p1,t1,s1] = ind2sub(size(MapProp_SiteLocB(:,targetsOnlyInProp,:)),linPB);
    linPBi = sub2ind(size(MapProp_SiteLoc),p1,targetsOnlyInProp(t1),s1,2*ones(length(t1),1));
    linPBo = sub2ind(size(MapQ_SiteLoc),p1,targetsOnlyInProp(t1),s1,2*ones(length(t1),1));
    MapQ_SiteLoc(linPBo) = MapProp_SiteLoc(linPBi);
    %MapQ_SiteLoc(:,targetsOnlyInCurr,1:numSitesCurr,:) = MapCurr_SiteLoc(:,targetsOnlyInCurr,1:numSitesCurr,:);  
    MapCurr_SiteLocA = squeeze(MapCurr_SiteLoc(:,:,:,1));
    linCA = find(MapCurr_SiteLocA(:,targetsOnlyInCurr,:));
    [p1,t1,s1] = ind2sub(size(MapCurr_SiteLocA(:,targetsOnlyInCurr,:)),linCA);
    linCAi = sub2ind(size(MapCurr_SiteLoc),p1,targetsOnlyInCurr(t1),s1,1*ones(length(t1),1));
    linCAo = sub2ind(size(MapQ_SiteLoc),p1,targetsOnlyInCurr(t1),s1,1*ones(length(t1),1));
    MapQ_SiteLoc(linCAo) = MapCurr_SiteLoc(linCAi);
    MapCurr_SiteLocB = squeeze(MapCurr_SiteLoc(:,:,:,2));
    linCB = find(MapCurr_SiteLocB(:,targetsOnlyInCurr,:));
    [p1,t1,s1] = ind2sub(size(MapCurr_SiteLocB(:,targetsOnlyInCurr,:)),linCB);
    linCBi = sub2ind(size(MapCurr_SiteLoc),p1,targetsOnlyInCurr(t1),s1,2*ones(length(t1),1));
    linCBo = sub2ind(size(MapQ_SiteLoc),p1,targetsOnlyInCurr(t1),s1,2*ones(length(t1),1));
    MapQ_SiteLoc(linCBo) = MapCurr_SiteLoc(linCBi);
    for v = 1:length(targetsOnlyInProp)
        Num_of_Molecule_Sites(targetsOnlyInProp(v)) = find(full(sum(squeeze(MapQ_DPS(:,targetsOnlyInProp(v),:)),1))==0,1)-1;
    end
    for v = 1:length(targetsOnlyInCurr)
        Num_of_Molecule_Sites(targetsOnlyInCurr(v)) = find(full(sum(squeeze(MapQ_DPS(:,targetsOnlyInCurr(v),:)),1))==0,1)-1;
    end 
    %KbMod
    linCP1 = find(MapProp_Vars{1}(:,targetsOnlyInProp,1:size(MapProp_Vars{1},3),1:size(MapProp_Vars{1},4)));
    [p1,t1,s1,m1] = ind2sub(size(MapProp_Vars{1}(:,targetsOnlyInProp,1:size(MapProp_Vars{1},3),1:size(MapProp_Vars{1},4))),linCP1);
    linCP1i = sub2ind(size(MapProp_Vars{1}),p1,targetsOnlyInProp(t1),s1,m1);
    linCP1o = sub2ind(size(MapQ_Vars{1}),p1,targetsOnlyInProp(t1),s1,m1);
    MapQ_Vars{1}(linCP1o) = MapProp_Vars{1}(linCP1i);
    linCV1 = find(MapCurr_Vars{1}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{1},3),1:size(MapCurr_Vars{1},4)));
    [p1,t1,s1,m1] = ind2sub(size(MapCurr_Vars{1}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{1},3),1:size(MapCurr_Vars{1},4))),linCV1);
    linCV1i = sub2ind(size(MapCurr_Vars{1}),p1,targetsOnlyInCurr(t1),s1,m1);
    linCV1o = sub2ind(size(MapQ_Vars{1}),p1,targetsOnlyInCurr(t1),s1,m1);
    MapQ_Vars{1}(linCV1o) = MapCurr_Vars{1}(linCV1i); 
    %TmMod
    linCP2 = find(MapProp_Vars{2}(:,targetsOnlyInProp,1:size(MapProp_Vars{2},3),1:size(MapProp_Vars{2},4)));
    [p2,t2,s2,m2] = ind2sub(size(MapProp_Vars{2}(:,targetsOnlyInProp,1:size(MapProp_Vars{2},3),1:size(MapProp_Vars{2},4))),linCP2);
    linCP2i = sub2ind(size(MapProp_Vars{2}),p2,targetsOnlyInProp(t2),s2,m2);
    linCP2o = sub2ind(size(MapQ_Vars{2}),p2,targetsOnlyInProp(t2),s2,m2);
    MapQ_Vars{2}(linCP2o) = MapProp_Vars{2}(linCP2i); 
    linCV2 = find(MapCurr_Vars{2}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{2},3),1:size(MapCurr_Vars{2},4)));
    [p2,t2,s2,m2] = ind2sub(size(MapCurr_Vars{2}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{2},3),1:size(MapCurr_Vars{2},4))),linCV2);
    linCV2i = sub2ind(size(MapCurr_Vars{2}),p2,targetsOnlyInCurr(t2),s2,m2);
    linCV2o = sub2ind(size(MapQ_Vars{2}),p2,targetsOnlyInCurr(t2),s2,m2);
    MapQ_Vars{2}(linCV2o) = MapCurr_Vars{2}(linCV2i);      
    %dHeqMod,dSeqMod,dCpMod
    linCP3A = find(MapProp_Vars{3}(:,targetsOnlyInProp,1:size(MapProp_Vars{3},3),1:size(MapProp_Vars{3},4)));
    linCP3B = find(MapProp_Vars{4}(:,targetsOnlyInProp,1:size(MapProp_Vars{4},3),1:size(MapProp_Vars{4},4)));
    linCP3C = find(MapProp_Vars{5}(:,targetsOnlyInProp,1:size(MapProp_Vars{5},3),1:size(MapProp_Vars{5},4)));
    linCP3 = unique([linCP3A linCP3B linCP3C]);
    [p3,t3,s3,m3] = ind2sub(size(MapProp_Vars{3}(:,targetsOnlyInProp,1:size(MapProp_Vars{3},3),1:size(MapProp_Vars{3},4))),linCP3);
    linCP3i = sub2ind(size(MapProp_Vars{3}),p3,targetsOnlyInProp(t3),s3,m3);
    linCP3o = sub2ind(size(MapQ_Vars{3}),p3,targetsOnlyInProp(t3),s3,m3);
    MapQ_Vars{3}(linCP3o) = MapProp_Vars{3}(linCP3i);
    MapQ_Vars{4}(linCP3o) = MapProp_Vars{4}(linCP3i);
    MapQ_Vars{5}(linCP3o) = MapProp_Vars{5}(linCP3i);
    linCV3A = find(MapCurr_Vars{3}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{3},3),1:size(MapCurr_Vars{3},4)));
    linCV3B = find(MapCurr_Vars{4}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{4},3),1:size(MapCurr_Vars{4},4)));
    linCV3C = find(MapCurr_Vars{5}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{5},3),1:size(MapCurr_Vars{5},4)));
    linCV3 = unique([linCV3A linCV3B linCV3C]);
    [p3,t3,s3,m3] = ind2sub(size(MapCurr_Vars{3}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{3},3),1:size(MapCurr_Vars{3},4))),linCV3);
    linCV3i = sub2ind(size(MapCurr_Vars{3}),p3,targetsOnlyInCurr(t3),s3,m3);
    linCV3o = sub2ind(size(MapQ_Vars{3}),p3,targetsOnlyInCurr(t3),s3,m3);
    MapQ_Vars{3}(linCV3o) = MapCurr_Vars{3}(linCV3i); 
    MapQ_Vars{4}(linCV3o) = MapCurr_Vars{4}(linCV3i); 
    MapQ_Vars{5}(linCV3o) = MapCurr_Vars{5}(linCV3i);
    %dHfMod,dSfMod,dHrMod,dSrMod
    linCP4A = find(MapProp_Vars{6}(:,targetsOnlyInProp,1:size(MapProp_Vars{6},3),1:size(MapProp_Vars{6},4)));
    linCP4B = find(MapProp_Vars{7}(:,targetsOnlyInProp,1:size(MapProp_Vars{7},3),1:size(MapProp_Vars{7},4)));
    linCP4C = find(MapProp_Vars{8}(:,targetsOnlyInProp,1:size(MapProp_Vars{8},3),1:size(MapProp_Vars{8},4)));
    linCP4D = find(MapProp_Vars{9}(:,targetsOnlyInProp,1:size(MapProp_Vars{9},3),1:size(MapProp_Vars{9},4)));
    linCP4 = unique([linCP4A linCP4B linCP4C linCP4D]);
    [p4,t4,s4,m4] = ind2sub(size(MapProp_Vars{6}(:,targetsOnlyInProp,1:size(MapProp_Vars{6},3),1:size(MapProp_Vars{6},4))),linCP4);
    linCP4i = sub2ind(size(MapProp_Vars{6}),p4,targetsOnlyInProp(t4),s4,m4);
    linCP4o = sub2ind(size(MapQ_Vars{6}),p4,targetsOnlyInProp(t4),s4,m4);
    MapQ_Vars{6}(linCP4o) = MapProp_Vars{6}(linCP4i);
    MapQ_Vars{7}(linCP4o) = MapProp_Vars{7}(linCP4i);
    MapQ_Vars{8}(linCP4o) = MapProp_Vars{8}(linCP4i);
    MapQ_Vars{9}(linCP4o) = MapProp_Vars{9}(linCP4i);
    linCV4A = find(MapCurr_Vars{6}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{6},3),1:size(MapCurr_Vars{6},4)));
    linCV4B = find(MapCurr_Vars{7}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{7},3),1:size(MapCurr_Vars{7},4)));
    linCV4C = find(MapCurr_Vars{8}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{8},3),1:size(MapCurr_Vars{8},4)));
    linCV4D = find(MapCurr_Vars{9}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{9},3),1:size(MapCurr_Vars{9},4)));
    linCV4 = unique([linCV4A linCV4B linCV4C linCV4D]);
    [p4,t4,s4,m4] = ind2sub(size(MapCurr_Vars{6}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{6},3),1:size(MapCurr_Vars{6},4))),linCV4);
    linCV4i = sub2ind(size(MapCurr_Vars{6}),p4,targetsOnlyInCurr(t4),s4,m4);
    linCV4o = sub2ind(size(MapQ_Vars{6}),p4,targetsOnlyInCurr(t4),s4,m4);
    MapQ_Vars{6}(linCV4o) = MapCurr_Vars{6}(linCV4i); 
    MapQ_Vars{7}(linCV4o) = MapCurr_Vars{7}(linCV4i); 
    MapQ_Vars{8}(linCV4o) = MapCurr_Vars{8}(linCV4i); 
    MapQ_Vars{9}(linCV4o) = MapCurr_Vars{9}(linCV4i);    
      %MapQ_Vars{2}(1,targetsOnlyInProp,1:size(MapProp_Vars{2},2)) = MapProp_Vars{2}(targetsOnlyInProp,1:size(MapProp_Vars{2},2));
%     linCP2 = find(MapProp_Vars{2}(targetsOnlyInProp,1:size(MapProp_Vars{2},2),1:size(MapProp_Vars{2},3)));
%     [t2,s2,m2] = ind2sub(size(MapProp_Vars{2}(targetsOnlyInProp,1:size(MapProp_Vars{2},2),1:size(MapProp_Vars{m2},3))),linCP2);
%     linCP2i = sub2ind(size(MapProp_Vars{2}),targetsOnlyInProp(t2),s2,m2);
%     linCP2o = sub2ind(size(MapQ_Vars{1}),1*ones(length(t2),1),targetsOnlyInProp(t2),s2,m2);
%     MapQ_Vars{2}(linCP2o) = MapProp_Vars{2}(linCP2i); 
%     if (map_idx==2)
%         %MapQ_Vars{2}(1,targetsOnlyInCurr,1:size(MapCurr_Vars{2},2)) = MapCurr_Vars{2}(targetsOnlyInCurr,1:size(MapCurr_Vars{2},2)); 
%         linCV2 = find(MapCurr_Vars{2}(targetsOnlyInCurr,1:size(MapCurr_Vars{2},2),1:size(MapCurr_Vars{2},3)));
%         [t2,s2,m2] = ind2sub(size(MapCurr_Vars{2}(targetsOnlyInCurr,1:size(MapCurr_Vars{2},2),1:size(MapCurr_Vars{2},3))),linCV2);
%         linCV2i = sub2ind(size(MapCurr_Vars{2}),targetsOnlyInCurr(t2),s2,m2);
%         linCV2o = sub2ind(size(MapQ_Vars{2}),1*ones(length(t2),1),targetsOnlyInCurr(t2),s2,m2);
%         MapQ_Vars{2}(linCV2o) = MapCurr_Vars{2}(linCV2i);  
%     else
%         linCV2 = find(MapCurr_Vars{2}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{2},3),1:size(MapCurr_Vars{2},4)));
%         [p2,t2,s2,m2] = ind2sub(size(MapCurr_Vars{2}(:,targetsOnlyInCurr,1:size(MapCurr_Vars{2},3),1:size(MapCurr_Vars{2},4))),linCV2);
%         linCV2i = sub2ind(size(MapCurr_Vars{2}),p2,targetsOnlyInCurr(t2),s2,m2);
%         linCV2o = sub2ind(size(MapQ_Vars{2}),p2,targetsOnlyInCurr(t2),s2,m2);
%         MapQ_Vars{2}(linCV2o) = MapCurr_Vars{2}(linCV2i);
%     end
    MapCurr_DPS = MapQ_DPS;
    MapCurr_SiteLoc = MapQ_SiteLoc;
    MapCurr_Vars = MapQ_Vars;
end
MapF_DPS = MapQ_DPS;
MapF_SiteLoc = MapQ_SiteLoc;
MapF_Vars = MapQ_Vars;
% MapF_Vars{2} = squeeze(max(MapQ_Vars{2},[],1));
end

% COMBINEBINDINGMAPS_PARFOR Combines two probe binding maps and their site locations for multiple targets
%   [Map3_DPS, Map3_SiteLoc] = COMBINEBINDINGMAPS_PARFOR(Map1_DPS, Map2_DPS, Map1_SiteLoc, Map2_SiteLoc) 
%   combines the binding maps in Map1_DPS and Map2_DPS and their corresponding site locations in Map1_SiteLoc and 
%   Map2_SiteLoc to generate Map3_DPS and Map3_SiteLoc for multiple targets. 



%plot plate
close all
outGroup2 = [1 4 8 9 11];
outGroup21 = [1 1 8 9 11];
outGroup = [1 4 8 9 11];
outGroup1 = [1 1 8 9 11];
outGroup3 = [1 4 8 9 11];
outGroup31 = [1 1 8 9 11];
Col = {[0 1 0],[0 0 0.8],[0 0 0.8],[0 0 0.8],[0 0 0.8],[0.4 0 0.6],[0.5 0 0.5],[1 0 1],[0 1 1],[0.9290 0.6940 0.1250],[0.9290 0.6940 0.1250]};opacity = 0.75;
theta = linspace(0,2*pi,100);r=0.5;
nPlates = 2;
T = readtable('Probes for pilot study final_ARF4.xlsm');
method = {'TS','SL','OS','PS','MF'};
TableStart = find(~cellfun(@isempty,T.PlateName_));
Probes_x0 = T.Position;
Probes_y0 = cellfun(@(x) find(ismember(flip(rows),x)),T.Plate);
Probes_Name = T.Sequence_5_To3___;
for nP = 1:length(TableStart)
figure(nP);
ax = gca;
rows = {'A','B','C','D','E','F','G','H'};
for xi = 1:12
    for yi = 1:8
    x0 = xi-0.5;y0 = yi-0.5;
    x = x0+r*cos(theta);
    y = y0+r*sin(theta);
    plot(x,y,'Color','k','LineWidth',0.5);hold on
    end
end
ax.XTick=[];
ax.YTick=[];
for xi = 1:12
   text(xi-0.5,8,num2str(xi),'vert','bottom','horiz','center','Color',[0 0 0],'FontSize',xFon,'FontWeight','Bold')             
end
for yi = 8:-1:1
   text(-0.5,yi-0.5,rows{9-yi},'vert','middle','horiz','center','Color',[0 0 0],'FontSize',xFon,'FontWeight','Bold')             
end
str = T.PlateName_(TableStart(nP));
text(0,9,replace(str,'_',' '),'vert','bottom','horiz','left','Color',[0 0 0],'FontSize',xFon,'FontWeight','Bold')            
ax.XLim = [0 12];
ax.YLim = [0 9];
    for v = 1:length(method)    
       text(12/8*v,-0.5,method{v},'vert','middle','horiz','center','Color',Col{outGroup2(v)},'FontSize',xFon,'FontWeight','Bold')       
    end
    text(12,-0.5,'Filtered Out','vert','middle','horiz','right','Color',[0.5 0.5 0.5],'FontSize',xFon,'FontWeight','Bold')        
end
for v = 1:length(method)
   Soft = find(contains(Probes_Name,method{v}));
   SoftT = Soft(contains(Probes_Name(Soft),'('));
   SoftK = setdiff(Soft,SoftT);
   for kk = 1:length(SoftK)
        figure(find(double(SoftK(kk)>=TableStart)==1,1,'last'));
        xi = Probes_x0(SoftK(kk));
        yi = Probes_y0(SoftK(kk));
        x0 = xi-0.5;y0 = yi-0.5;
        x = x0+r*cos(theta);
        y = y0+r*sin(theta);
        patch(x,y,Col{outGroup(v)},'FaceAlpha',1);
        text(x0,y0,extractAfter(Probes_Name(SoftK(kk)),strcat(method{v},'-')),'vert','middle','horiz','center','Color',[0 0 0],'FontSize',xFon,'FontWeight','Bold')                
   end
   for kk = 1:length(SoftT)
        figure(find(double(SoftT(kk)>=TableStart)==1,1,'last'));
        xi = Probes_x0(SoftT(kk));
        yi = Probes_y0(SoftT(kk));
        x0 = xi-0.5;y0 = yi-0.5;
        x = x0+r*cos(theta);
        y = y0+r*sin(theta);
        patch(x,y,Col{outGroup(v)},'FaceAlpha',1);
        text(x0,y0,extractBetween(extractAfter(Probes_Name(SoftT(kk)),strcat(method{v},'-')),'(',')'),'vert','middle','horiz','center','Color',[0.5 0.5 0.5],'FontSize',xFon,'FontWeight','Bold')                  
   end
end
for nP = 1:length(TableStart)
figure(nP);
str = T.PlateName_(TableStart(nP));
saveas(gcf,strcat(replace(str{:},'/','_'),'_PlateLayout.png'))
end
%0 1
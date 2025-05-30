ui = 1;vi = 7;
illumination_corrected = 1;
int_level = 1;
local_over_cell_bkg = 1;
useFixedThreshold = 0;
get_figdata

% plot intensity & RNA count distributions Ver 5
% clear all;
% load('figData.mat');
sp = 1;
fs = 10;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'};
figure(20); clf;
for rep = 1:3;
    xbin = [0:20:1000];
    for soft = 2:5;
        subplot(6,4,sp);
        ydata = figData.signal_values{soft,rep};
        ydataRef = figData.signal_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef);
        pAllInt(soft,rep) = p;
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([300 700]); ylim([0 0.02]);
        xlabel('Intensity'); 
        if soft == 2;
            ylabel('Probability');
        elseif soft == 5;
            % Add replica label manually using 'text'
            ax = gca;
            yl = ylim(ax);
            xPos = ax.Position(1) + ax.Position(3) + 0.005; % slightly right of current axis
            yPos = ax.Position(2) + ax.Position(4)/4;      % vertically centered

            annotation('textbox', [xPos, yPos, 0.05, 0.05], ...
                       'String', ['Rep ' num2str(rep)], ...
                       'EdgeColor', 'none', ...
                       'HorizontalAlignment', 'left', ...
                       'VerticalAlignment', 'middle', ...
                       'FontSize', fs, ...
                       'FontWeight', 'bold');
        else;
        end;
        if rep == 1;
            title([softw{soft}]);
        else
        end;
        set(gca, 'FontSize', fs);
        sp = sp + 1;
    end;
end;
%legend('TrueProbes', softw{soft});

for rep = 1:3;
    xbin = [0:1:20];
    for soft = 2:5;
        subplot(6,4,sp);
        ydata = figData.spot_count_values{soft,rep};
        ydataRef = figData.spot_count_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef);
        pAllRNA(soft,rep) = p;
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([-1 20]); ylim([0 0.4]);
        xlabel('No. of RNA molecules');
        if soft == 2;
            ylabel('Probability');
        elseif soft == 5;
            % Add replica label manually using 'text'
            ax = gca;
            yl = ylim(ax);
            xPos = ax.Position(1) + ax.Position(3) + 0.005; % slightly right of current axis
            yPos = ax.Position(2) + ax.Position(4)/4;      % vertically centered

            annotation('textbox', [xPos, yPos, 0.05, 0.05], ...
                       'String', ['Rep ' num2str(rep)], ...
                       'EdgeColor', 'none', ...
                       'HorizontalAlignment', 'left', ...
                       'VerticalAlignment', 'middle', ...
                       'FontSize', fs, ...
                       'FontWeight', 'bold');

        else;
        end;
        %title(['p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        sp = sp + 1;
    end;
end;
%legend('TrueProbes', softw{soft});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all;
% load('figData.mat');

% plot intensity & RNA count distributions Ver 1
lw = 5;
figure(1); clf;
for rep = 1:3;
    subplot(2,3,rep);
    xbin = [0:20:1000];
    for soft = 1:5;
        ydata = figData.signal_values{soft,rep};
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs',LineWidth=lw); hold on;
    end;
    xlim([300 700]);
    ylim([0 0.03]);
    xlabel('Intensity'); ylabel('Probability');

    subplot(2,3,rep+3);
    xbin = [0:1:20];
    for soft = 1:5;
        ydata = figData.spot_count_values{soft,rep};
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs',LineWidth=lw); hold on;
    end;
    xlim([-2 20]);
    ylim([0 0.5]);
    xlabel('No. of RNA molecules'); ylabel('Probability');
end;
legend('TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH');


% plot intensity & RNA count distributions Ver 2
figure(2); clf;
sp = 1;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'}
for rep = 1:3;
    xbin = [0:20:1000];
    for soft = 1:5;
        subplot(3,5,sp);
        ydata = figData.signal_values{soft,rep};
        ydataRef = figData.signal_values{1,rep};
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        sp = sp + 1;
        xlim([300 700]);
        ylim([0 0.03]);
        title(['Rep: ' num2str(rep) ', ' softw{soft}]);
    end;
end;

    subplot(2,3,rep+3);
    xbin = [0:1:20];
    for soft = 1:5;
        ydata = figData.spot_count_values{soft,rep};
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs',LineWidth=5); hold on;
    end;
    xlim([0 25]);
    ylim([0 1]);
%end;
legend('TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH');


% plot intensity & RNA count distributions Ver 3
%clear all;
%load('figData.mat');
sp = 1;
fs = 20;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'}
for rep = 1:3;
    figure(3+rep); clf;
    xbin = [0:20:1000];
    for soft = 2:5;
        subplot(2,4,soft-1);
        ydata = figData.signal_values{soft,rep};
        ydataRef = figData.signal_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([300 700]); ylim([0 0.02]);
        xlabel('Intensity'); ylabel('Probability');
        title(['Rep: ' num2str(rep) ', ' softw{soft} ', p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;

    xbin = [0:1:20];
    for soft = 2:5;
        subplot(2,4,soft+4-1);
        ydata = figData.spot_count_values{soft,rep};
        ydataRef = figData.spot_count_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([0 20]); ylim([0 0.4]);
        xlabel('No. of RNA molecules'); ylabel('Probability');
        title(['p=' num2str(p)]);
        set(gca, 'FontSize', fs);
                legend(softw{soft},'TrueProbes');
    end;
end;

% plot intensity & RNA count distributions Ver 4
%clear all;
%load('figData.mat');
sp = 1;
fs = 20;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'};
ls = {'-','--','-.'};
numData = 3;
shadesBlue = [linspace(0.2, 0, numData)', linspace(0.2, 0, numData)', linspace(0.8, 1, numData)'];
shadesRed = [linspace(1, 0.5, numData)', linspace(0.2, 0, numData)', linspace(0.2, 0, numData)'];
figure(10); clf;
soft = 2;
for soft = 2:5;
    ydataAll = NaN(2000,3);
    ydataRefAll = ydataAll;
    subplot(2,4,soft-1);
    for rep = 1:3;
        xbin = [0:20:1000];
        ydata = figData.signal_values{soft,rep};
        ydataAll(1:length(ydata),rep) = ydata;
        ydataRef = figData.signal_values{1,rep};
        ydataRefAll(1:length(ydataRef),rep) = ydataRef;       
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5,'EdgeColor',shadesRed(rep,:)); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5,'EdgeColor',shadesBlue(rep,:)); hold on;
    end;
    ydataAll1 = ydataAll(:);
    ydataAll1 = ydataAll1(~isnan(ydataAll1));
    ydataRefAll1 = ydataRefAll(:);
    ydataRefAll1 = ydataRefAll1(~isnan(ydataRefAll1));
    [h,p] = kstest2(ydataAll1,ydataRefAll1);
    xlim([300 700]); ylim([0 0.02]);
    xlabel('Intensity'); ylabel('Probability');
    set(gca, 'FontSize', fs);
            legend(softw{soft},'TrueProbes');
    title(['Rep: ' num2str(rep) ', ' softw{soft} ', p=' num2str(p)]);

    ydataAll = NaN(2000,3);
    ydataRefAll = ydataAll;
    subplot(2,4,soft+4-1);
    for rep = 1:3;
        xbin = [0:1:20];
        ydata = figData.spot_count_values{soft,rep};
        ydataAll(1:length(ydata),rep) = ydata;
        ydataRef = figData.spot_count_values{1,rep};
        ydataRefAll(1:length(ydataRef),rep) = ydataRef;       
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5,'EdgeColor',shadesRed(rep,:)); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5,'EdgeColor',shadesBlue(rep,:)); hold on;
    end;
    ydataAll1 = ydataAll(:);
    ydataAll1 = ydataAll1(~isnan(ydataAll1));
    ydataRefAll1 = ydataRefAll(:);
    ydataRefAll1 = ydataRefAll1(~isnan(ydataRefAll1));
    [h,p] = kstest2(ydataAll1,ydataRefAll1);
    xlim([0 20]); ylim([0 0.4]);
    xlabel('No. of RNA molecules'); ylabel('Probability');
    title(['p=' num2str(p)]);
    set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
end;

% plot intensity & RNA count distributions Ver 3
%clear all;
%load('figData.mat');
sp = 1;
fs = 14;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'}
for rep = 1:3;
    figure(40+rep); clf;
    xbin = [0:20:1000];
    for soft = 2:5;
        subplot(2,4,soft-1);
        ydata = figData.backgd_values{soft,rep};
        ydataRef = figData.backgd_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([200 500]); %ylim([0 0.02]);
        xlabel('Background Intensity'); ylabel('Probability');
        title(['Rep: ' num2str(rep) ', ' softw{soft} ', p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;

    xbin = [0:20:1000]
    for soft = 2:5;
        subplot(2,4,soft+4-1);
        ydata = figData.signal_minus_backgd_values{soft,rep};
        ydataRef = figData.signal_minus_backgd_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([0 200]);% ylim([0 0.4]);
        xlabel('Signal Minus Background'); ylabel('Probability');
        title(['p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;
end;



sp = 1;
fs = 14;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH'}
for rep = 1:3;
    figure(60+rep); clf;
    xbin = [0:20:1000];
    for soft = 2:5;
        subplot(2,4,soft-1);
        ydata = figData.backgd_values{soft,rep};
        ydataRef = figData.backgd_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([200 500]); %ylim([0 0.02]);
        xlabel('Background Intensity'); ylabel('Probability');
        title(['Rep: ' num2str(rep) ', ' softw{soft} ', p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;

    xbin = [0:1:40]
    for soft = 2:5;
        subplot(2,4,soft+4-1);
        ydata = figData.SNR_values{soft,rep};
        ydataRef = figData.SNR_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([0 10]);% ylim([0 0.4]);
        xlabel('Signal-to-Noise Ratio (SNR)'); ylabel('Probability');
        title(['p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;
end;

sp = 1;
fs = 14;
softw = {'TrueProbes', 'Stellaris', 'Oligostan', 'PaintSHOP', 'MERFISH','No Probe'}
for rep = 1:3;
    figure(80+rep); clf;
    xbin = [0:20:1000];
    for soft = 2:6;
        subplot(2,5,soft-1);
        ydata = figData.cell_backgd_mean_values{soft,rep};
        ydataRef = figData.cell_backgd_mean_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([200 360]);% ylim([0 0.02]);
        xlabel('Cell Background Intensity'); ylabel('Probability');
        title(['Rep: ' num2str(rep) ', ' softw{soft} ', p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;

    xbin = [0:1:40]
    for soft = 2:6;
        subplot(2,5,soft+5-1);
        ydata = figData.cell_backgd_fano_values{soft,rep};
        ydataRef = figData.cell_backgd_fano_values{1,rep};
        [h,p] = kstest2(ydata,ydataRef)
        h = histogram(ydata,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        h1 = histogram(ydataRef,xbin,'Normalization','pdf','DisplayStyle','stairs','LineWidth',5); hold on;
        xlim([0 9]); %ylim([0 0.4]);
        xlabel('Cell Background Fano-Factor'); ylabel('Probability');
        title(['p=' num2str(p)]);
        set(gca, 'FontSize', fs);
        legend(softw{soft},'TrueProbes');
    end;
end;
min(cellfun(@(x) min(x), figData.signal_values),[],'all')
max(cellfun(@(x) max(x), figData.signal_values),[],'all')
min(cellfun(@(x) min(x), figData.backgd_values),[],'all')
max(cellfun(@(x) max(x), figData.backgd_values),[],'all')
min(cellfun(@(x) min(x), figData.signal_minus_backgd_values),[],'all')
max(cellfun(@(x) max(x), figData.signal_minus_backgd_values),[],'all')
min(cellfun(@(x) min(x), figData.SNR_values),[],'all')
max(cellfun(@(x) max(x), figData.SNR_values),[],'all')
min(cellfun(@(x) min(x), figData.spot_count_values),[],'all')
max(cellfun(@(x) max(x), figData.spot_count_values),[],'all')
min(cellfun(@(x) min(x), figData.cell_backgd_mean_values),[],'all')
max(cellfun(@(x) max(x), figData.cell_backgd_mean_values),[],'all')
min(cellfun(@(x) min(x), figData.cell_backgd_fano_values),[],'all')
max(cellfun(@(x) max(x), figData.cell_backgd_fano_values),[],'all')


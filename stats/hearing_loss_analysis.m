%%%
%   NeuroTin / 20231014 / hearing_loss_stats.m
%
%   Re-assess HL differences between the CBT and fMRI groups
%   as suggested (revision).
%%%
%
%   Nicolas Gninenko / nicolas.gninenko@gmail.com
%
%%%


% Load data
clear; clc;
%load('TinnitusAlloc_GnVdata_NG.mat');
%load('HL_TinnLevels.mat');
%load('CBT_fMRI_patients_cells.mat');
load('231011 matlab.mat'); % everything from workspace directly


%% Plot CPT-AMAs (average between Left and Right, in %)
% for both groups, in order to visualize that the CBT group
% has one outlier with nearly 60% of HL
close all;
tmp_f1 = figure;
plot(TinnitusAllocGnVdata_CBT{:,8},'LineWidth',1.2); hold on;
plot(TinnitusAllocGnVdata_fMRI{:,8},'LineWidth',1.2);
tmp_f1.CurrentAxes.XLim = [0 23]; tmp_f1.CurrentAxes.XTick = 1:22;
tmp_f1.CurrentAxes.YLim = [0 65]; tmp_f1.CurrentAxes.YTick = 0:10:60;
tmp_f1.CurrentAxes.FontSize = 10;
%tmp_f1.CurrentAxes.YGrid = 'on';
tmp_f1.CurrentAxes.FontName = 'Basis Grotesque Pro';
title('Average CPT-AMA (%) for CBT and fMRI groups at baseline','FontSize',16,...
    'FontName','Basis Grotesque Pro','FontWeight','Normal');
xlabel('Participants','FontSize',14,'FontName','Basis Grotesque Pro');
ylabel('%','FontSize',14,'FontName','Basis Grotesque Pro');
plot([1 22],[20 20],'black','Color',[0 0 0 0.25]);
legend([{['CBT (mean: ' num2str(mean(TinnitusAllocGnVdata_CBT{:,8})) ...
    ' ± ' num2str(std(TinnitusAllocGnVdata_CBT{:,8})) ' [SD])']},{...
    ['fMRI (mean: ' num2str(mean(TinnitusAllocGnVdata_fMRI{:,8})) ' ± ' ...
    num2str(std(TinnitusAllocGnVdata_fMRI{:,8})) ' [SD])']}],'Location','Best',...
    'Box','on','FontSize',16,'FontName','Basis Grotesque Pro');
text(4.75,40,{['CBT (median: ' num2str(median(TinnitusAllocGnVdata_CBT{:,8})) ' ± ' ...
    num2str(mad(TinnitusAllocGnVdata_CBT{:,8},1)) ' [MAD])'],...
    ['fMRI (median: ' num2str(median(TinnitusAllocGnVdata_fMRI{:,8})) ' ± ' ...
    num2str(mad(TinnitusAllocGnVdata_fMRI{:,8},1)) ' [MAD])']},...
    'FontSize',16,'FontName','Basis Grotesque Pro');


%% Plot CPT-AMA (average between Left and Right, in %)
% for both groups without the outlier in CBT

close all;
tmp_f2 = figure;
plot(TinnitusAllocGnVdata_CBT{[1:8 10:end],8},'LineWidth',1.2); hold on;
plot(TinnitusAllocGnVdata_fMRI{:,8},'LineWidth',1.2);
tmp_f2.CurrentAxes.XLim = [0 23]; tmp_f2.CurrentAxes.XTick = 1:22;
tmp_f2.CurrentAxes.YLim = [0 65]; tmp_f2.CurrentAxes.YTick = 0:10:60;
tmp_f2.CurrentAxes.FontSize = 10;
%tmp_f2.CurrentAxes.YGrid = 'on';
tmp_f2.CurrentAxes.FontName = 'Basis Grotesque Pro';
title('Average CPT-AMA (%) for CBT and fMRI groups','FontSize',16,...
    'FontName','Basis Grotesque Pro','FontWeight','Normal');
xlabel('Participants','FontSize',14,'FontName','Basis Grotesque Pro');
ylabel('%','FontSize',14,'FontName','Basis Grotesque Pro');
plot([1 22],[20 20],'black','Color',[0 0 0 0.25]);
legend([{['CBT (mean: ' num2str(mean(TinnitusAllocGnVdata_CBT{[1:8 10:end],8})) ...
    ' ± ' num2str(std(TinnitusAllocGnVdata_CBT{[1:8 10:end],8})) ' [SD])']},{...
    ['fMRI (mean: ' num2str(mean(TinnitusAllocGnVdata_fMRI{:,8})) ' ± ' ...
    num2str(std(TinnitusAllocGnVdata_fMRI{:,8})) ' [SD])']}],'Location','Best',...
    'Box','on','FontSize',16,'FontName','Basis Grotesque Pro');
text(4.75,40,{['CBT (median: ' num2str(median(TinnitusAllocGnVdata_CBT{[1:8 10:end],8})) ' ± ' ...
    num2str(mad(TinnitusAllocGnVdata_CBT{[1:8 10:end],8},1)) ' [MAD])'],...
    ['fMRI (median: ' num2str(median(TinnitusAllocGnVdata_fMRI{:,8})) ' ± ' ...
    num2str(mad(TinnitusAllocGnVdata_fMRI{:,8},1)) ' [MAD])']},...
    'FontSize',16,'FontName','Basis Grotesque Pro');


%% Bilateral HL at Tinn. Frequency for Left and Right ears

tmp_f3 = figure;
subplot(211);
plot(sort(CBT_HearingLoss_TinnLevels{:,2}),'LineWidth',1.2); hold on;
plot(sort(fMRI_HearingLoss_TinnLevels{:,2}),'LineWidth',1.2); hold off;
tmp_f3.CurrentAxes.XLim = [0 23]; tmp_f3.CurrentAxes.XTick = 1:22;
tmp_f3.CurrentAxes.YLim = [0 110]; tmp_f3.CurrentAxes.YTick = 0:20:120;
tmp_f3.CurrentAxes.FontSize = 10; tmp_f3.CurrentAxes.YGrid = 'on';
tmp_f3.CurrentAxes.FontName = 'Basis Grotesque Pro';
title('Hearing loss values at tinnitus frequency at baseline','FontSize',16,...
    'FontName','Basis Grotesque Pro','FontWeight','Normal');
ylabel('dB HL','FontSize',14,'FontName','Basis Grotesque Pro');
legend({'CBT','fMRI'},'Location','Best',...
    'Box','on','FontSize',14,'FontName','Basis Grotesque Pro');
tmp_pval = num2str(ranksum(fMRI_HearingLoss_TinnLevels{:,2},CBT_HearingLoss_TinnLevels{:,2}));
text(3,80,['P = ' tmp_pval ' (Mann-Whitney U test)'],...
    'Interpreter','tex','FontSize',12,'FontName','Basis Grotesque Pro');
text(13,20,'Left ear','FontSize',12,'FontName','Basis Grotesque Pro');
subplot(212);
plot(sort(CBT_HearingLoss_TinnLevels{:,4}),'LineWidth',1.2); hold on;
plot(sort(fMRI_HearingLoss_TinnLevels{:,4}),'LineWidth',1.2); hold off;
tmp_f3.CurrentAxes.XLim = [0 23]; tmp_f3.CurrentAxes.XTick = 1:22;
tmp_f3.CurrentAxes.YLim = [0 110]; tmp_f3.CurrentAxes.YTick = 0:20:120;
tmp_f3.CurrentAxes.FontSize = 10; tmp_f3.CurrentAxes.YGrid = 'on';
tmp_f3.CurrentAxes.FontName = 'Basis Grotesque Pro';
%title('Sorted hearing loss [dB HL] values at tinnitus frequency (right ear)','FontSize',14,...
%    'FontName','Basis Grotesque Pro','FontWeight','Normal');
ylabel('dB HL','FontSize',14,'FontName','Basis Grotesque Pro');
xlabel('Participants (sorted)','FontSize',14,'FontName','Basis Grotesque Pro');
tmp_pval = num2str(ranksum(fMRI_HearingLoss_TinnLevels{:,4},CBT_HearingLoss_TinnLevels{:,4}));
text(3,80,['P = ' tmp_pval ' (Mann-Whitney U test)'],...
    'Interpreter','tex','FontSize',12,'FontName','Basis Grotesque Pro');
text(13,20,'Right ear','FontSize',12,'FontName','Basis Grotesque Pro');


%% Bilateral Tinn. Perception at the given Frequency(ies) of Hearing Loss


tmp_f4 = figure;
%hold off;
subplot(211);
plot(sort(CBT_HearingLoss_TinnLevels{:,3}-CBT_HearingLoss_TinnLevels{:,2}),'LineWidth',1.2); hold on;
plot(sort(fMRI_HearingLoss_TinnLevels{:,3}-fMRI_HearingLoss_TinnLevels{:,2}),'LineWidth',1.2); hold off;
tmp_f4.CurrentAxes.XLim = [0 23]; tmp_f4.CurrentAxes.XTick = 1:22;
tmp_f4.CurrentAxes.YLim = [0 40]; tmp_f4.CurrentAxes.YTick = 5:5:35;
tmp_f4.CurrentAxes.FontSize = 10; tmp_f4.CurrentAxes.YGrid = 'on';
tmp_f4.CurrentAxes.FontName = 'Basis Grotesque Pro';
title('Tinnitus loudness w.r.t. hearing loss at tinn. frequency at baseline','FontSize',16,...
    'FontName','Basis Grotesque Pro','FontWeight','Normal');
ylabel('dB SL','FontSize',14,'FontName','Basis Grotesque Pro');
tmp_pval = num2str(ranksum(fMRI_HearingLoss_TinnLevels{:,3}-fMRI_HearingLoss_TinnLevels{:,2},...
    CBT_HearingLoss_TinnLevels{:,3}-CBT_HearingLoss_TinnLevels{:,2}));
text(2,25,['P = ' tmp_pval ' (Mann-Whitney U test)'],...
    'Interpreter','tex','FontSize',12,'FontName','Basis Grotesque Pro');
text(17,10,'Left ear','FontSize',12,'FontName','Basis Grotesque Pro');
subplot(212);
plot(sort(CBT_HearingLoss_TinnLevels{:,5}-CBT_HearingLoss_TinnLevels{:,4}),'LineWidth',1.2); hold on;
plot(sort(fMRI_HearingLoss_TinnLevels{:,5}-fMRI_HearingLoss_TinnLevels{:,4}),'LineWidth',1.2); hold off;
tmp_f4.CurrentAxes.XLim = [0 23]; tmp_f4.CurrentAxes.XTick = 1:22;
tmp_f4.CurrentAxes.YLim = [0 40]; tmp_f4.CurrentAxes.YTick = 5:5:35;
tmp_f4.CurrentAxes.FontSize = 10; tmp_f4.CurrentAxes.YGrid = 'on';
tmp_f4.CurrentAxes.FontName = 'Basis Grotesque Pro';
%title('Sorted hearing loss [dB HL] values at tinnitus frequency (right ear)','FontSize',14,...
%    'FontName','Basis Grotesque Pro','FontWeight','Normal');
ylabel('dB SL','FontSize',14,'FontName','Basis Grotesque Pro');
xlabel('Participants (sorted)','FontSize',14,'FontName','Basis Grotesque Pro');
tmp_pval = num2str(ranksum(fMRI_HearingLoss_TinnLevels{:,5}-fMRI_HearingLoss_TinnLevels{:,4},...
    CBT_HearingLoss_TinnLevels{:,5}-CBT_HearingLoss_TinnLevels{:,4}));
%tmp_pval = tmp_pval(2:4);
text(2,25,['P = ' tmp_pval ' (Mann-Whitney U test)'],...
    'Interpreter','tex','FontSize',12,'FontName','Basis Grotesque Pro');
text(17,10,'Right ear','FontSize',12,'FontName','Basis Grotesque Pro');
legend({'CBT','fMRI'},'Location','NorthEast',...
    'Box','on','FontSize',14,'FontName','Basis Grotesque Pro');
%export_fig[...]


%% Levene's quad. test to assess homoscedasticity

[p_,stats] = vartestn([TinnitusAllocGnVdata_fMRI{:,9};...
TinnitusAllocGnVdata_fMRI{:,10};...
TinnitusAllocGnVdata_CBT{:,9};...
TinnitusAllocGnVdata_CBT{:,10}],...
[ones(2*length(TinnitusAllocGnVdata_fMRI{:,9}),1);...
2*ones(2*length(TinnitusAllocGnVdata_CBT{:,9}),1)],'TestType','LeveneQuadratic');
% Bartlett's
[h_,p_] = vartest2([TinnitusAllocGnVdata_fMRI{:,9};TinnitusAllocGnVdata_fMRI{:,10}],...
    [TinnitusAllocGnVdata_CBT{:,9};TinnitusAllocGnVdata_CBT{:,10}]);


%% Non-parametric comparisons for PTA

[h_,p_] = ranksum(TinnitusAllocGnVdata_CBT{:,9},TinnitusAllocGnVdata_fMRI{:,9});   % .0516 -> .052
[h_,p_] = ranksum(TinnitusAllocGnVdata_CBT{:,10},TinnitusAllocGnVdata_fMRI{:,10}); % .0661 -> .066


%% Non-parametric and parametric comparisons for PTA without the outlier (participant 9 of the CBT group)

[h_,p_] = ranksum(TinnitusAllocGnVdata_CBT{[1:8 10:22],9},TinnitusAllocGnVdata_fMRI{:,9});
[h_,p_] = ranksum(TinnitusAllocGnVdata_CBT{[1:8 10:22],10},TinnitusAllocGnVdata_fMRI{:,10});

[h_,p_] = ttest2(TinnitusAllocGnVdata_CBT{[1:8 10:22],9},TinnitusAllocGnVdata_fMRI{:,9});
[h_,p_] = ttest2(TinnitusAllocGnVdata_CBT{[1:8 10:22],10},TinnitusAllocGnVdata_fMRI{:,10});




%% Principal Component Analysis
% Note: This file saves at several points, hit any key to resume

%% Input and Load Parameters

% Path to emiily_functions
addpath(fullfile('~','Documents','MATLAB','emily_functions','projection'))
path2ems_functions = fullfile('~','Documents','MATLAB','emily_functions');
addpath(path2ems_functions)

% Path to input files
addpath(fullfile('~','Dropbox','Data','DC_GEO_Data','jake'))
inputMat = '/Users/caz3so/Dropbox/Data/Miraldi/rotation_methods_data/DC_GEO_Data/jake/kk21_geneExprMeta.mat';
[path2data,matName,ext] = fileparts(inputMat);

% Output Figure Folder
figFolder = fullfile(path2data,'PCA');
mkdir(figFolder)

% Load Functions and Parameters
load cbcolors.mat
load linecolors.mat
colors = cbcolors;
colors(8,:) = mean([cbcolors(8,:);zeros(1,3)]);
colors(6,:) = [];
colors(end,:) = [];
colors = [linecolors(9,:); colors; linecolors(18:19,:)];
load(fullfile(path2data,matName))
dataset = matName;
printDname = [strrep(dataset,'_',' ')];

%% PCA Settings

% Normalization
dataScale = 'log2';
    % 'zscore' -- mean-center, variance normalize **NOTE: log2 won't work with CLR-normalized data due to negatives
    % 'log2' -- add a pseudocount, divide by the mean, log2-normalize
    % '' -- no data transformation will be done, raw input quantification
    % 'sqzscore' -- mean-center, sqrt variance normalize

% Plot Settings
metaCatOfInt = 'Group'; %'Sample	Group	Treatment	Sorting
sampleLabelCat = 'Treatment'; % Label Names
loadOn = 1;  % Loadings?  1 = yes, 0 = no
totWeights = 50; % Number of genes plotted in loadings plot
labOn = 1; % Scores dots on loading plot?  1 = yes, 0 = no

%% The Magic

[~, obs] = size(ncounts);
currconds = ones(obs,1);
inLimGroup = 1:obs;

metaInd = find(ismember(sampMetaCats,metaCatOfInt));
currcond_nums2 = currconds;
groupLabels = 'Treatment';
if metaInd
    currGroups = {sampMetaData{:,metaInd}};
    uniGroups = unique(currGroups);
    totGroups = length(uniGroups);
    for gp = 1:totGroups
        currGroup = uniGroups{gp};
        currcond_nums2(ismember(currGroups,currGroup)) = gp;
    end
    groupLabels = replace_underscore_w_spaces( uniGroups');       
end

%% Normalize
[vars, obs] = size(ncounts);

alg = '';

if length(find(ismember({dataScale},'log2')))
    dataMean = mean(ncounts+1,2);
    dataScaled = log2((ncounts+1)./repmat(dataMean,1,obs));
    typename = 'log2(FC)';
    disp(typename)
elseif length(find(ismember({dataScale},'zscore')))
    dataScaled = zscore(ncounts')';
    typename = 'z-score';
    disp(typename)
elseif length(find(ismember({dataScale},'sqzscore')))
    dataScaled = sqrtZscore(ncounts')';
    typename = 'sqrt(z-score)';
    disp(typename)    
else
    dataScaled = ncounts;
    typename = 'Raw';
    disp(typename)
end

currdata = dataScaled;

pca_type = 'Classic';
alpha = .95;
out_cutoff = .95;
plots = [1 1 1 1];
numweights = 10000;

frac_measure_minimum = size(currdata,1);
plotweights = 1:min(numweights,size(currdata,1));

group_nums = unique(currcond_nums2);
totgroups = length(group_nums);

[coefs, scores, latent, tsquared, var_exp] = princomp(currdata','econ');
size(scores)
pcs = 4;

%% Scores Plot with Ellipse -- PCs 1& 2

%Figure Set up
figure(1), clf
out_cutoff = .95;
axis auto

%Plot Points- Marker, Color, and Size selections
if and(labOn,length(groupLabels))
    hold on
    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),'o',...
            'Color',colors(group_ind,:),'MarkerSize',16, 'MarkerFaceColor', colors(group_ind,:))
    end
    legend(groupLabels,'FontSize',16)
end

% for group_ind = 1:totgroups
%     group = group_nums(group_ind);
 %    text(scores(currcond_nums2==group, 1),scores(currcond_nums2==group,2),...
 %        strvcat(printSnames{currcond_nums2==group}),...
 %        'Color',colors(group_ind,:),'FontWeight','Bold',...
 %        'FontSize',14)
 %end

plot(1.1*min(scores(:,1)),1.2*min(scores(:,2)),'w')
hold on
plot(1.2*max(scores(:,1)),1.1*max(scores(:,2)),'w')    
ax = axis();
plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5],'LineWidth',3)
plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5],'LineWidth',3)
axis manual

stddev2 = sqrt((std(scores(:,2))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
stddev1 = sqrt((std(scores(:,1))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));

%Axis Label Options
title([printDname ' ' typename ' Scores (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
    'FontSize',22)
xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
x = [-stddev1:stddev1/50:stddev1];
yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
yl = -yu;
plot(x,yu,'Color',[.5 .5 .5],'LineWidth',3)
plot(x,yl,'Color',[.5 .5 .5],'LineWidth',3)
hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')
set(gca,'XTick',[],'YTick',[],'LineWidth',2), box on

if find(ismember({metaCatOfInt},'SIZE'))
    colormap(jet), colorbar('Location','EastOutside','YTick',[0 1],'YTickLabel',...
        {'LowCov','High Cov'},'FontSize',14)
end
%% Loadings Plot PC1/2
if loadOn
    norms = sqrt(sum(coefs(:,1:2).^2,2));
    [normsort, sortinds] = sort(norms,'descend');
    toPlot = sortinds(1:totWeights);
    load_scale = min([ax(2)/max(coefs(coefs(:,1)>0,1))...
        ax(1)/min(coefs(coefs(:,1)<0,1))...
        ax(4)/max(coefs(coefs(:,2)>0,2)),...
        ax(3)/min(coefs(coefs(:,2)<0,2))]);

    figure(2), clf
    hold on     

    if and(labOn,length(groupLabels))
        hold on
        for group_ind = 1:totgroups
            group = group_nums(group_ind);
            plot(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),'.',...
                'Color',colors(group_ind,:),'MarkerSize',20)
        end
        legend(groupLabels,'FontSize',16)
    end

    set(gca,'XTick',[],'YTick',[],'LineWidth',2), box on
    plot(load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,1)]',...
        load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,2)]',...
        'y','LineWidth',4)
    hold on
    ax = axis();
    plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5])
    plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5])
    plot(load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,1)]',...
        load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,2)]',...
        'y','LineWidth',4)
    
    if and(labOn,length(groupLabels))
        hold on
        for group_ind = 1:totgroups
            group = group_nums(group_ind);
            plot(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),'.',...
                'Color',colors(group_ind,:),'MarkerSize',20)
        end
    end
    
    text(load_scale*coefs(toPlot,1)',...
        load_scale*coefs(toPlot,2)',cellstr(strvcat(genesc{toPlot})),'FontWeight','Bold')
    title([strrep(dataset,'_',' ') ' ' typename ' Loadings (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
        'FontSize',24)
    xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
        'FontSize',20)
    ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
        'FontSize',20)
    if find(ismember({metaCatOfInt},'SIZE'))
        colormap(jet), colorbar('Location','EastOutside','YTick',[0 1],'YTickLabel',...
            {'LowCov','High Cov'},'FontSize',14)
    end
end 

%% Scores Plot with Ellipse -- PCs 3 & 4
figure(3), clf
out_cutoff = .95;
axis auto

if and(labOn,length(groupLabels))
    hold on
    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),'o',...
            'Color',colors(group_ind,:),'MarkerSize',16, 'MarkerFaceColor', colors(group_ind,:))
    end
    legend(groupLabels,'FontSize',16)
end

plot(1.1*min(scores(:,3)),1.2*min(scores(:,4)),'w')
hold on
plot(1.2*max(scores(:,3)),1.1*max(scores(:,4)),'w')        
ax = axis();
plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5],'LineWidth',3)
plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5],'LineWidth',3)
axis manual
stddev2 = sqrt((std(scores(:,4))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
stddev1 = sqrt((std(scores(:,3))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
title([printDname ' ' typename ' Scores (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
    'FontSize',24)
xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
x = [-stddev1:stddev1/50:stddev1];
yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
yl = -yu;
plot(x,yu,'Color',[.5 .5 .5],'LineWidth',3)
plot(x,yl,'Color',[.5 .5 .5],'LineWidth',3)
hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')

set(gca,'XTick',[],'YTick',[],'LineWidth',2), box on

if find(ismember({metaCatOfInt},'SIZE'))
    colormap(jet), colorbar('Location','EastOutside','YTick',[0 1],'YTickLabel',...
        {'LowCov','High Cov'},'FontSize',16)
end
%% Loadings Plot PC3/4
if loadOn
    norms = sqrt(sum(coefs(:,3:4).^2,2));
    [normsort, sortinds] = sort(norms,'descend');
    toPlot = sortinds(1:totWeights);

    load_scale = min([ax(2)/max(coefs(coefs(:,3)>0,1))...
        ax(1)/min(coefs(coefs(:,3)<0,1))...
        ax(4)/max(coefs(coefs(:,4)>0,2)),...
        ax(3)/min(coefs(coefs(:,4)<0,2))]);     

    figure(3), clf

    if and(labOn,length(groupLabels))
        hold on
        for group_ind = 1:totgroups
            group = group_nums(group_ind);
            plot(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),'.',...
                'Color',colors(group_ind,:),'MarkerSize',16)
        end
        legend(groupLabels,'FontSize',16)
    end

    plot(load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,3)]',...
        load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,4)]',...
        'y','LineWidth',4)

    hold on
    ax = axis();
    plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5])
    plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5])
    plot(load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,3)]',...
        load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,4)]',...
        'y','LineWidth',4)
    
    if and(labOn,length(groupLabels))
        hold on
        for group_ind = 1:totgroups
            group = group_nums(group_ind);
            plot(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),'.',...
                'Color',colors(group_ind,:),'MarkerSize',20)
        end
    end
    text(load_scale*coefs(toPlot,3)',...
        load_scale*coefs(toPlot,4)',cellstr(strvcat(genesc{toPlot})),'FontWeight','Bold')
    title([strrep(dataset,'_',' ') ' ' typename ' Loadings (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
        'FontSize',24)
    xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
        'FontSize',20)
    ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
        'FontSize',20)
    set(gca,'XTick',[],'YTick',[],'LineWidth',2), box on
    if find(ismember({metaCatOfInt},'SIZE'))
        colormap(jet), colorbar('Location','EastOutside','YTick',[0 1],'YTickLabel',...
            {'LowCov','High Cov'},'FontSize',14)
    end
end
    
%% PLOT 3 axis at once!

figure(200),clf
if and(labOn,length(groupLabels))
    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot3(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
            scores(currcond_nums2==group,3),'o',...
            'Color',colors(group_ind,:),'MarkerSize',16, 'MarkerFaceColor', colors(group_ind,:))
            hold on
    end   
    legend(groupLabels,'FontSize',16)
end

plot3(1.1*min(scores(:,1)),1.2*min(scores(:,2)),1.2*min(scores(:,3)),'w')
hold on
plot3(1.2*max(scores(:,1)),1.1*max(scores(:,2)),1.2*max(scores(:,3)),'w')    
plot3(scores(:,1),scores(:,2),scores(:,3),'w.')
for samp = 1: length(scores(:,1))
    plot3([0 scores(samp,1)],[0,scores(samp,2)],[0,scores(samp,3)],'y-')
end

hold on

if find(ismember({metaCatOfInt},'SIZE'))
    colormap(jet), colorbar('Location','EastOutside','YTick',[0 1],'YTickLabel',...
        {'LowCov','High Cov'},'FontSize',14)
end

grid on

ax = axis();
plot3([ax(1) ax(2)],[0 0],[0 0],'Color',[.5 .5 .5],'LineWidth',3)
plot3([0 0],[ax(3) ax(4)],[0 0],'Color',[.5 .5 .5],'LineWidth',3)
plot3([0 0],[0 0],[ax(5) ax(6)],'Color',[.5 .5 .5],'LineWidth',3)
xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
zlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
currfig = [figinf '_pcs3D_scores'];
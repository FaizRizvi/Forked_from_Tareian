%% pca_anal_geneExp_meta_limGenes_limGroups_ivTh_1710
%% can limit to a set of genes 
%% Also, this file helps visualize:
% 1. size factors (from DESeq2 normalization) / coverage
% 2. metadata
%% Note: this file s at several points so that you can adjust the 
% figure size on the monitor before it saves to the adjusted dimensions.
% You can also move labels to improve figure legibility.
% Hit any key to resume, once you're satisfied with adjusted dimensions.

%% begin parameters

addpath(fullfile('~','erm','MATLAB','Genesets'))
addpath(fullfile('~','erm','MATLAB','emily_functions','projection'))
path2ems_functions = fullfile('~','erm','MATLAB','emily_functions');

% datatag = 'Tr1_Th17_biasLM'; % will be used when saving figure outputs

%% end parameters

% matName = [datatag];
% dataset = strrep(matName,'_noBatch.mat','');

inputMat = '/Users/emiraldi/erm/tmpRfolder/lymph36broad/lymph36broad_VSD.mat';

% dataScale = 'log2';%'log2'; % 'log2'; % '', 'zscore' are options
dataScale = ''; % Set the data scaling options; current possibilities
    % are **NOTE: log2 won't work with CLR-normalized data due to negatives
    % 'zscore' -- mean-center, variance normalize
    % 'log2' -- add a pseudocount, divide by the mean, log2-normalize
    % '' -- no data transformation will be done, raw input quantification
    %       will be used
    % 'sqzscore' -- mean-center, sqrt variance normalize

loadOn = 1;  % plot loadings?  1 = yes, 0 = no

totWeights = 50; % number of genes plotted in loadings

labOn = 0; % add scores dots to loading plot?  1 = yes, 0 = no


for limGroupInd = 3 % which groups of samples to limit --> 3 = ivTh

printOpt = 0; % whether to output a table of counts corresponding to the groups and sets of genes of interest 

% How to color your samples, based on sampMetaCats
metaCatOfInt = 'ivThCondition';%'broadCond';%'broadCond';%'seqExp';%SIZE'; % 
%% choices include: 'condition, broadCond, seqExp, SIZE, conditionMinusBatch, ivThCondition
% 'SIZE' loads size factors

sampleLabelCat = 'conditionMinusBatch'; % what to use as names for labeling samples

% If you want to remove samples
removeSamples = {...'DKO_Ileum_ILC3_CCR6+';
	'Th0_DMSO_Donor_1_day_5_SL5251';
	'Th17_DMSO_Donor_1_day_5_SL5252';
	'Th0_DMSO_Donor_2_day_5_SL5254';
	'Th17_DMSO_Donor_2_day_5_SL5255';
	'Th17_C16_Donor_1_day_5_SL5253';
	'Th17_C16_Donor_2_day_5_SL5256'};

% If you decided to limit genes, files and names, can be left blank
limitGeneLists = {'';
    '/Users/emiraldi/erm/tmpRfolder/Th17_48h_RNAseq/GeneSets/union_Th17genesMM10_Th17vTh0_FC0p58_FDR10.txt';
    '/Users/emiraldi/erm/tmpRfolder/Th17_48h_RNAseq/GeneSets/Th17Paper_genesOnly.txt';
    '/Users/emiraldi/erm/tmpRfolder/Th17_48h_RNAseq/GeneSets/Th17paper_48h_up_FC0p58_FDR25.txt';
    '/Users/emiraldi/erm/tmpRfolder/Th17_48h_RNAseq/GeneSets/Th17genes_48h_up_FC0p58_FDR25.txt'};
limitGeneInfs = {'';
    'currentTargs';
    '17pap';
    '17union';
    '17me'};

limGeneInd = 2; % the list of genes (e.g., Th17) to use

limitGenes = limitGeneLists{limGeneInd};
limitGeneInf = limitGeneInfs{limGeneInd};

% If you'd like to limit by cell type
limGroups = {{'Bcell';
    'ILC_LI';
    'ILC_SI';
    'Th_LI';
    'Th_MLN';
    'Th_SC';
    'Th_SI';
    'Th_Spl';
    'evThymocyte';
    'ivTh'};
    {'ivTh';
    'Th_LI';
    'Th_MLN';
    'Th_SC';
    'Th_SI';
    'Th_Spl'};
    {'ivTh'}};
limGroupInfs = {'';'allTh';'ivTh'};
limGroupMeta = 'broadCond'; % <-- in the metadata this is the name of the metadata category, sampMetaCats

limGroup = limGroups{limGroupInd};
limGroupInf = limGroupInfs{limGroupInd};



[path2data,matName,ext] = fileparts(inputMat);
if length(limGroup) 
    figFolder = fullfile(path2data,['PCA_anal_' limGroupInf]);
else
    figFolder = fullfile(path2data,['PCA_anal' limGroupInf]);
end

addpath(path2ems_functions) % add path to functions
load cbcolors.mat
load linecolors.mat
colors = cbcolors;
colors(8,:) = mean([cbcolors(8,:);zeros(1,3)]);
colors(6,:) = [];
colors(end,:) = [];
colors = [linecolors(9,:); colors; linecolors(18:19,:)];

mkdir(figFolder)

%% PCA analysis
load(fullfile(path2data,matName))
dataset = matName
printDname = [strrep(dataset,'_',' ')];

[vars obs] = size(ncounts);

currconds = ones(obs,1);

%% LIMIT to sample conditions
if length(limGroupInf)
    limMetaInd = find(ismember(sampMetaCats,limGroupMeta));
    allMetaGroups = {sampMetaData{:,limMetaInd}};
    inLimGroup = [];
    totLimGroups = length(limGroup);
    for lg = 1:totLimGroups
        inLimGroup = [inLimGroup;
            find(ismember(allMetaGroups,limGroup{lg}))'];
    end
else
    inLimGroup = 1:obs;
end

% note: you could color samples according to group membership, which would
% be encoded in currcond_nums2 (a unique number for each group that would
% be ordered according to the labels in "printSnames").  I didn't want to
% hardcode this in, so there are no colorings for this example file.
metaInd = find(ismember(sampMetaCats,metaCatOfInt));
currcond_nums2 = currconds;
groupLabels = '';
if metaInd
    currGroups = {sampMetaData{:,metaInd}};
    uniGroups = unique(currGroups);
    totGroups = length(uniGroups);
    for gp = 1:totGroups
        currGroup = uniGroups{gp};
        currcond_nums2(find(ismember(currGroups,currGroup))) = gp;
    end
    groupLabels = replace_underscore_w_spaces( uniGroups');
%% fun stuff if you have numerical metadata, e.g., sequencing depth
    elseif find(ismember({metaCatOfInt},'SIZE'))
        [vals inds] = sort(sizeFactors);
        currcond_nums2 = inds;
        colors = colormap((jet(length(currcond_nums2))));        
end

%% limit genes, if necessary
if length(limitGenes) > 0
    fin = fopen(limitGenes,'r');
    C = textscan(fin,'%s');
    [vals, keepInds] = intersect(lower(genesc),lower(C{1}));
    genesc = cellstr(strvcat(genesc{keepInds}));
    ncounts = ncounts(keepInds,:);
    disp(['Genes Limited to ' limitGenes])
    disp([num2str(length(keepInds)) ' genes.'])
    [tt,limFil,ext] = fileparts(limitGenes);
    bitInf = {limFil};  
    figinf = fullfile(figFolder,bitInf{1},[matName '_' limitGeneInf '_' metaCatOfInt '_' dataScale ]);
    mkdir(fullfile(figFolder,bitInf{1}))
else
    figinf = fullfile(figFolder,[matName '_' metaCatOfInt '_' dataScale ]);
end

%% limit samples to the groups of interest (e.g., ivTh) & pick correct sample name
sampleLabInd = find(ismember(sampMetaCats,sampleLabelCat));
sampLabels = {sampMetaData{inLimGroup,sampleLabInd}}';
conditionsc = cellstr(strvcat(conditionsc{inLimGroup}));
ncounts = ncounts(:,inLimGroup);
currcond_nums2 = currcond_nums2(inLimGroup);
[vars obs] = size(ncounts);

%% remove any samples (e.g., outliers specified above)
tot2remove = length(removeSamples);
if tot2remove
    disp(['Tot to remove = ' num2str(tot2remove)])
    removeInd = [];
    for samp = 1:tot2remove
        ind = find(ismember(conditionsc,removeSamples{samp}));
        removeInd = [removeInd; ind];
    end
    disp(['REMOVED ' num2str(length(find(removeInd))) ':'])
    strvcat(conditionsc{removeInd})

    [vals, keep] = setdiff(1:obs,removeInd);
    conditionsc = cellstr(strvcat(conditionsc{keep}));
    ncounts = ncounts(:,keep);
    currcond_nums2 = currcond_nums2(keep);
    sampLabels = cellstr(strvcat(sampLabels{keep}));
end 

% remove some group labels if there's no longer anything in those
% categories
if length(find(ismember({'SIZE'},metaCatOfInt))) == 0
    uniGroups = unique(currcond_nums2);
    totGroups = length(uniGroups);
    groupLabels = cellstr(strvcat(groupLabels{uniGroups}));
end
% gompers

sampleNames = sampLabels;

printSnames = replace_underscore_w_spaces(cellstr(sampleNames));

%% print data to table, if desired
if printOpt
    printFile = [strrep(figinf,['_' dataScale],'') '_' limGroupInf '.txt'];
    fid = fopen(printFile,'w');
    fprintf(fid,['\t' strjoin(conditionsc,'\t') '\n']);
    for gene = 1:vars
        fprintf(fid,[genesc{gene} '\t' strjoin(cellstr(num2str(ncounts(gene,:)')),'\t') '\n']);
    end
    fclose(fid);
    disp(printFile)
    size(ncounts)
end

end

% gompers

%% Normalize
[vars obs] = size(ncounts);

alg = '';


%% Apply data normalization
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
figure(1), clf
out_cutoff = .95;
% subplot(1,10,2:5)
%% PC plane 1 X 2
% get axis set
axis auto
if and(labOn,length(groupLabels))
    hold on
    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),'.',...
            'Color',colors(group_ind,:),'MarkerSize',20)
    end
    legend(groupLabels,'FontSize',16)
end
plot(1.1*min(scores(:,1)),1.2*min(scores(:,2)),'w')
hold on
plot(1.2*max(scores(:,1)),1.1*max(scores(:,2)),'w')    
% plot(scores(:,1),scores(:,2),'k','MarkerSize',12)
ax = axis();
plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5],'LineWidth',3)
plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5],'LineWidth',3)
axis manual
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
        strvcat(printSnames{currcond_nums2==group}),...
        'Color',colors(group_ind,:),'FontWeight','Bold',...
        'FontSize',12)
end

stddev2 = sqrt((std(scores(:,2))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
stddev1 = sqrt((std(scores(:,1))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
title([printDname ' ' typename ' Scores (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
    'FontSize',24)
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


pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');


% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs12_scores.pdf'])
currfig = [figinf '_pcs12_scores'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])

 
%% Loadings
% subplot(1,10,[7:10])
% totWeights = 75;
if loadOn
    norms = sqrt(sum(coefs(:,1:2).^2,2));
    [normsort sortinds] = sort(norms,'descend');
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
    % plot(load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,1)]',...
    %     load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,2)]',...
    %     'y','LineWidth',4)
    plot(load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,1)]',...
        load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,2)]',...
        'y','LineWidth',4)

    hold on
    ax = axis();
    plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5])
    plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5])
    % plot(load_scale*[zeros(size(coefs(:,1))) coefs(:,1)]',...
    %     load_scale*[zeros(size(coefs(:,1))) coefs(:,2)]',...
    %     'y','LineWidth',4)
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
%         legend(groupLabels,'FontSize',16)
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
    
    pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

    % fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
    % print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
    currfig = [figinf '_pcs12_loadings'];
    save2pdf([currfig '.pdf'],gcf,150)
    saveas(gcf,currfig,'fig')
    disp([currfig '.pdf'])
end 

%% Scores Plot with Ellipse -- PCs 3 & 4
figure(3), clf
out_cutoff = .95;
% subplot(1,10,2:5)
%% PC plane 3, 4
% get axis set
axis auto

if and(labOn,length(groupLabels))
    hold on
    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),'.',...
            'Color',colors(group_ind,:),'MarkerSize',20)
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
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    text(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),...
        replace_underscore_w_spaces(cellstr(strvcat(printSnames{currcond_nums2==group}))),...
        'Color',colors(group_ind,:),'FontWeight','Bold',...
        'FontSize',12)
%     text(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),...
%         'o','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
%     text(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),...
%         'x','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
end
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
        {'LowCov','High Cov'},'FontSize',14)
end

pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
currfig = [figinf '_pcs34_scores'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])

%% PC Plane 3/4 Loadings
if loadOn
    norms = sqrt(sum(coefs(:,3:4).^2,2));
    [normsort sortinds] = sort(norms,'descend');
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
                'Color',colors(group_ind,:),'MarkerSize',20)
        end
        legend(groupLabels,'FontSize',16)
    end

    % plot(load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,1)]',...
    %     load_scale*[zeros(size(coefs(toPlot,1))) coefs(toPlot,2)]',...
    %     'y','LineWidth',4)
    plot(load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,3)]',...
        load_scale*[zeros(size(coefs(toPlot,3))) coefs(toPlot,4)]',...
        'y','LineWidth',4)

    hold on
    ax = axis();
    plot([ax(1) ax(2)],[0 0],'Color',[.5 .5 .5])
    plot([0 0],[ax(3) ax(4)],'Color',[.5 .5 .5])
    % plot(load_scale*[zeros(size(coefs(:,1))) coefs(:,1)]',...
    %     load_scale*[zeros(size(coefs(:,1))) coefs(:,2)]',...
    %     'y','LineWidth',4)
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
%         legend(groupLabels,'FontSize',16)
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
    
    pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

    % fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
    % print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
    currfig = [figinf '_pcs34_loadings'];
    save2pdf([currfig '.pdf'],gcf,150)
    saveas(gcf,currfig,'fig')
    disp([currfig '.pdf'])
end
    
%% PLOT 3 axis at once!

figure(200),clf
if and(labOn,length(groupLabels))

    for group_ind = 1:totgroups
        group = group_nums(group_ind);
        plot3(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
            scores(currcond_nums2==group,3),'.',...
            'Color',colors(group_ind,:),'MarkerSize',20)
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
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    
    text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
        scores(currcond_nums2==group,3),...
        strvcat(printSnames{currcond_nums2==group}),...
        'Color',colors(group_ind,:),'FontWeight','Bold',...
        'FontSize',14)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'o','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'x','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
end

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
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf + .fig'])

% %% Pareto Diagram -- Variance Explained
% 
% T = scores;
% P = coefs;
% denom = sum(diag(dataset*dataset'));
% var_exp = zeros(pcs,1);
% cum_var_exp = zeros(pcs,1);
% for pc = 1:pcs
%     t = T(:,pc);
%     p = P(:,pc);
%     var_exp(pc) = sum(diag((t*p')*(t*p')'))/denom;
%     if pc > 1
%         cum_var_exp(pc) = var_exp(pc) + cum_var_exp(pc-1);
%     else
%         cum_var_exp(pc) = var_exp(pc);
%     end
% end
% 
% figure(10), clf
% bar(100*var_exp)
% hold on
% plot(100*cum_var_exp,'bo')
% plot(100*cum_var_exp,'b')
% set(gca,'XTick',1:pcs,'XTickLabel',num2str([1:pcs]'),...
%     'FontWeight','Bold')
% axis tight
% xlabel('Principle Components','FontWeight','Bold')
% ylabel('% Varience Explained','FontWeight','Bold')
% title('Variance Explained','FontWeight','Bold')


% 
% %% Loadings
% % subplot(1,10,[7:10])
% norms = sqrt(sum(coefs.^2,2));
% [normsort sortinds] = sort(norms,'descend');
% plotweights = sortinds(plotweights);
% 
% figure(2), clf
% load_scale = min([ax(2)/max(coefs(coefs(:,1)>0,1))...
%     ax(1)/min(coefs(coefs(:,1)<0,1))...
%     ax(4)/max(coefs(coefs(:,2)>0,2)),...
%     ax(3)/min(coefs(coefs(:,2)<0,2))]);     
% % plot(load_scale*[zeros(size(coefs(plotweights,1))) coefs(plotweights,1)]',...
% %     load_scale*[zeros(size(coefs(plotweights,1))) coefs(plotweights,2)]',...
% %     '.','MarkerSize',10,'Color',[.5 .5 .5])
% plot(load_scale*coefs(plotweights,1)',...
%     load_scale*coefs(plotweights,2)',...
%     '.','MarkerSize',10,'Color',[.5 .5 .5])
% 
% 
% hold on
% 
% 
% colorinds = [4 3 1];
% 
% for se = 1:totsets
%     currset = weightgroups(:,se);
%     currset = currset(plotweights);
%     disp(weightsetnames{se+1})
%     length(find(currset))
%     plot(load_scale*currset'.*coefs(plotweights,1)',...
%         load_scale*currset'.*coefs(plotweights,2)',...
%         '.','MarkerSize',20,'Color',colors(colorinds(se),:))
% end
% 
% set(gca,'FontSize',20,'FontWeight','Bold')
% legend(weightsetnames)
% 
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% axis tight
% ax = axis();
% plot([ax(1) ax(2)],[0 0],'k')
% plot([0 0],[ax(3) ax(4)],'k')
% 
% title([printDname ' ' typename ' Loadings (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
%     'FontSize',24)
% xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% 
% 
% % for 
% 
% % 
% % text(load_scale*coefs(plotweights,1)',...
% %     load_scale*coefs(plotweights,2)',cellstr(strvcat(currlabels{plotweights})),'FontWeight','Bold')
% % title([printDname ' ' typename ' Loadings (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
% %     'FontSize',24)
% % xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
% %     'FontSize',20)
% % ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
% %     'FontSize',20)
% 
% pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');
% 
% % fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% % print('-painters','-dpdf','-r1800',[figinf '_pcs12_loadings.pdf'])
% save2pdf([figinf '_pcs12_loadings.pdf'],gcf,300)
% saveas(gcf,[figinf '_pcs12_loadings'],'fig')
% 
% %% Loadings
% % % % subplot(1,10,[7:10])
% figure(4), clf
% load_scale = min([ax(2)/max(coefs(coefs(:,3)>0,3))...
%     ax(1)/min(coefs(coefs(:,3)<0,3))...
%     ax(4)/max(coefs(coefs(:,4)>0,4)),...
%     ax(3)/min(coefs(coefs(:,4)<0,4))]);     
% plot(load_scale*coefs(plotweights,3)',...
%     load_scale*coefs(plotweights,4)',...
%     '.','MarkerSize',10,'Color',[.5 .5 .5])
% hold on
% 
% for se = 1:totsets
%     currset = weightgroups(:,se);
%     currset = currset(plotweights);
%     disp(weightsetnames{se+1})
%     length(find(currset))
%     plot(load_scale*currset'.*coefs(plotweights,3)',...
%         load_scale*currset'.*coefs(plotweights,4)',...
%         '.','MarkerSize',20,'Color',colors(colorinds(se),:))
% end
% 
% set(gca,'FontSize',20,'FontWeight','Bold')
% legend(weightsetnames)
% 
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% axis tight
% ax = axis();
% plot([ax(1) ax(2)],[0 0],'k')
% plot([0 0],[ax(3) ax(4)],'k')
% 
% title([printDname ' ' typename ' Loadings (' alg ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
%     'FontSize',24)
% xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% 
% pause, disp('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');
% 
% % fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% % print('-painters','-dpdf','-r1800',[figinf '_pcs34_loadings.pdf'])
% save2pdf([figinf '_pcs34_loadings.pdf'],gcf,300)
% saveas(gcf,[figinf '_pcs34_loadings'],'fig')
% 
% %% PC plane 3 X 4
% % % disp('entered')
% % figure
% % subplot(3,6,13)
% % axis auto
% % plot(min(scores(:,3)),min(scores(:,4)),'y')
% % hold on
% % plot(max(scores(:,3)),max(scores(:,4)),'y')
% % axis tight
% % 
% % ax = axis();
% % plot([ax(1) ax(2)]',[0 0]','r')
% % plot([0 0]',[ax(3) ax(4)]','r')
% % axis manual
% % for group_ind = 1:totgroups
% %     group = group_nums(group_ind);
% %     text(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),...
% %         cellstr(num2str(currlivernums{currcond_nums2==group})),'Color',...
% %         colors(group_ind,:),'FontWeight','Bold')
% % end
% % plot(scores(~outr.flag.sd,3),scores(~outr.flag.sd,4),'m*')
% % plot(scores(~outr.flag.od,3),scores(~outr.flag.od,4),'g*')
% % plot(scores(and(~outr.flag.od,~outr.flag.sd),1),...
% %     scores(and(~outr.flag.od,~outr.flag.sd),2),'r*')
% % stddev2 = sqrt((std(scores(:,4))^2)*finv(out_cutoff,pcs,obs-pcs)*2*obs/(obs-pcs));
% % stddev1 = sqrt((std(scores(:,3))^2)*finv(out_cutoff,pcs,obs-pcs)*2*obs/(obs-pcs));
% % x = [-stddev1:stddev1/50:stddev1];
% % yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
% % yl = -yu;
% % plot(x,yu)
% % plot(x,yl)
% % hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
% % hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')
% % title([pca_type 'Scores'],'FontWeight','Bold')
% % xlabel('PC 3','FontWeight','Bold')
% % ylabel('PC 4','FontWeight','Bold')
% % end
% 
% %% Scores Plane Distance and Orthogonal Distance to the Plane
% % if plots(2)
% %     subplot(3,6,14)
% %     plot(outr.sd,outr.od,'b*')
% %     axis tight
% %     ax=axis();
% %     hold on,plot([outr.cutoff.sd outr.cutoff.sd],[ax(3) ax(4)],'r')
% %     hold on,plot([ax(1) ax(2)],[outr.cutoff.od outr.cutoff.od],'r')
% %     plot(outr.sd(~outr.flag.sd),outr.od(~outr.flag.sd),'m*')
% %     plot(outr.sd(~outr.flag.od),outr.od(~outr.flag.od),'g*')
% %     plot(outr.sd(and(~outr.flag.od,~outr.flag.sd)),...
% %         outr.od(and(~outr.flag.od,~outr.flag.sd)),'r*')
% %     text(outr.sd(~outr.flag.sd)',outr.od(~outr.flag.sd)',...
% %         strvcat(currlivernumsc(~outr.flag.sd)),'FontWeight','Bold')%,'Color','m')
% %     text(outr.sd(~outr.flag.od),outr.od(~outr.flag.od),...
% %         strvcat(currlivernumsc(~outr.flag.od)),'FontWeight','Bold')%,'g*')
% %     text(outr.sd(and(~outr.flag.od,~outr.flag.sd)),...
% %         outr.od(and(~outr.flag.od,~outr.flag.sd)),...
% %         strvcat(currlivernumsc(and(~outr.flag.od,~outr.flag.sd))),...
% %         'FontWeight','Bold')%,'r*')
% %     title('Position','FontWeight','Bold')
% %     xlabel('Distance within  PC Plane','FontWeight','Bold')
% %     ylabel('Orthogonal Distance to PC Plane','FontWeight','Bold')
% % % end
% % 
% % 
% % 
% % % if pcs > 1
% % %     figure
% % %         plot(scores(:,1),scores(:,2),'b*')
% % %         %gname(currlivernums)
% % %         hold on
% % % 
% % %         plot(scores(~outr.flag.sd,1),scores(~outr.flag.sd,2),'m*')
% % %         plot(scores(~outr.flag.od,1),scores(~outr.flag.od,2),'g*')
% % %         plot(scores(and(~outr.flag.od,~outr.flag.sd),1),...
% % %             scores(and(~outr.flag.od,~outr.flag.sd),2),'r*')
% % %         stddev2 = sqrt((std(scores(:,2))^2)*finv(out_cutoff,2,obs-2)*2*obs/(obs-2));
% % %         stddev1 = sqrt((std(scores(:,1))^2)*finv(out_cutoff,2,obs-2)*2*obs/(obs-2));
% % % 
% % %         title(pca_type,'FontWeight','Bold')
% % %         xlabel('PC 1','FontWeight','Bold')
% % %         ylabel('PC 2','FontWeight','Bold')
% % %         x = [-stddev1:stddev1/50:stddev1];
% % %         yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
% % %         yl = -yu;
% % %         plot(x,yu)
% % %         plot(x,yl)
% % %         axis tight
% % %         ax = axis();
% % %         hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % %         hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % %         hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
% % %         hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')
% % %         text(scores(:,1)',...
% % %             scores(:,2)',currlivernums)
% % % %         text(max(stddev1,stddev2)*coefs(:,1)',...
% % % %             max(stddev1,stddev2)*coefs(:,2)',currlabels)
% % % end
% % % 
% % % if pcs > 1
% % %     figure
% % %         plot(scores(:,3),scores(:,4),'b*')
% % %         %gname(currlivernums)
% % %         hold on
% % % % 
% % % %         plot(scores(~outr.flag.sd,1),scores(~outr.flag.sd,2),'m*')
% % % %         plot(scores(~outr.flag.od,1),scores(~outr.flag.od,2),'g*')
% % % %         plot(scores(and(~outr.flag.od,~outr.flag.sd),1),...
% % % %             scores(and(~outr.flag.od,~outr.flag.sd),2),'r*')
% % % %         stddev2 = sqrt((std(scores(:,2))^2)*finv(out_cutoff,2,obs-2)*2*obs/(obs-2));
% % % %         stddev1 = sqrt((std(scores(:,1))^2)*finv(out_cutoff,2,obs-2)*2*obs/(obs-2));
% % % 
% % %         title(pca_type,'FontWeight','Bold')
% % %         xlabel('PC 3','FontWeight','Bold')
% % %         ylabel('PC 4','FontWeight','Bold')
% % %         %zlabel('PC 3','FontWeight','Bold')
% % % %         x = [-stddev1:stddev1/50:stddev1];
% % % %         yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
% % % %         yl = -yu;
% % % %         plot(x,yu)
% % % %         plot(x,yl)
% % %         axis tight
% % % %         ax = axis();
% % % %         hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % % %         hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
% % % %         hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
% % % %         hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')
% % %         text(scores(:,3)',...
% % %             scores(:,4)',scores(:,3),currlivernums)
% % % %         text(max(stddev1,stddev2)*coefs(:,1)',...
% % % %             max(stddev1,stddev2)*coefs(:,2)',currlabels)
% 
% 



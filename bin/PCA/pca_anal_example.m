%% pca_anal_example
%% Goals:
% 1. load the text file
% 2. filter data
% 3. normalize the data
% 4. visualize PCA scores and loading plots
% When you go through this code, enter each command one line at a time to
% see what it does.  Or for lines that apply multiple functions run parts
% of the line at a time to see what they do.

%% These provide links to some of my custom functions; please change to 
% so that the path corresponds to the location of my custom functions on
% your computer
addpath(fullfile('~','Documents','MATLAB','emily_functions','projection'))
addpath(fullfile('~','Documents','MATLAB','emily_functions'))
addpath(fullfile('~','Dropbox','Data','DC_GEO_Data','jake'))
%% Begin parameters
% Note: you'll need to change the path to the gene expression file to
% wherever it is on your computer
inputGeneExpressionFile = 'KK221_expression.txt';
outputFolder = 'kk21_pca';
dataset = 'kk21_pca';  % add a dataset name, this will be used in figure titles and output files

loadingOpt = 1;     % Plot loadings?  1 --> yes, 0 --> no
totWeights = 75;    % What number of top genes would you like to see?
    % Looking at all ~10k results in a blob.
dataScale = 'log2'; % Set the data scaling options; current possibilities

geneExpressionCutoff = 1; % filter genes, if they do not reach this minimum 
visualizeGeneExpressionCut = 1; % visualize histogram of sample maximum 
fontSize = 12; % font size for figures    

%% End parameters, Begin analysis
    
load linecolors.mat     % My custom color scheme matrix, saved as a .mat in emily_functions
mkdir(outputFolder)        % Make the outputFolder, in case it doesn't already exist
figinf = fullfile(outputFolder,[dataset '_' dataScale ]);  % Make a base name for output figures

%% Load text file
% First see how many conditions / columns are in the file by reading the
% first line with tgetl
fid = fopen(inputGeneExpressionFile,'r');
tline = fgetl(fid);
fclose(fid);
sampleNames = cellstr(strvcat(strsplit(tline,'\t')));    % get sample names
totSamps = length(sampleNames(2:end)); 
% Now that we know # of samples, reopen and get the rest
fid = fopen(inputGeneExpressionFile,'r');
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',1);
fclose(fid);
genesc = C{1};          % get gene names
ncounts = [C{2:end}];   % get gene expression values

% create figure text variables with underscores removed to make figure text look better
printDname = [strrep(dataset,'_',' ')];
printSnames = replace_underscore_w_spaces(cellstr(sampleNames));

[vars obs] = size(ncounts);

%% Apply gene expression cutoff
disp(['Applying gene expression cutoff: ' num2str(geneExpressionCutoff) '.'])
geneMaxes = max(ncounts,[],2);
% visualize maximum gene expression and cutoff if you want
if visualizeGeneExpressionCut
    figure(100),clf
    hist(geneMaxes,1000) 
    title([printDname ' maximum gene expression values'],'FontSize',fontSize+2)
    xlabel('Maximum expression level per gene','FontSize',fontSize)
    ylabel('Counts','FontSize',fontSize)
    axis tight
    % plot cutoff
    hold on    
    ax = axis();
    plot(geneExpressionCutoff*[1 1],ax(3:4),'r','LineWidth',2)
    legend({'Counts','Cutoff'},'Location','Best')
    set(gca,'FontSize',fontSize)
end
highEnoughInds = find(geneMaxes>=geneExpressionCutoff);
genesc = cellstr(strvcat(genesc{highEnoughInds}));
ncounts = ncounts(highEnoughInds,:);

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
else
    dataScaled = ncounts;
    typename = 'Raw';
    disp(typename)
end

plots = [1 1 1 1];
numweights = 10000;

% make sure that totWeights isn't bigger than number of genes
totGenes = size(dataScaled,1);
plotweights = 1:min(totWeights,size(dataScaled,1));


% note: you could color samples according to group membership, which would
% be encoded in currcond_nums2 (a unique number for each group that would
% be ordered according to the labels in "printSnames").  I didn't want to
% hardcode this in, so there are no colorings for this example file.
currcond_nums2 = ones(obs,1);
group_nums = unique(currcond_nums2);
totgroups = length(group_nums);

[coefs, scores, latent, tsquared, var_exp] = princomp(dataScaled','econ');
size(scores)
pcs = 4;

%% Scores Plot with Ellipse -- PCs 1& 2
figure(1), clf
out_cutoff = .95;
% subplot(1,10,2:5)
%% PC plane 1 X 2
% get axis set
axis auto
plot(1.1*min(scores(:,1)),1.2*min(scores(:,2)),'w')
hold on
plot(1.2*max(scores(:,1)),1.1*max(scores(:,2)),'w')    
% plot(scores(:,1),scores(:,2),'k','MarkerSize',12)
ax = axis();
plot([ax(1) ax(2)],[0 0],'Color',linecolors(8,:),'LineWidth',3)
plot([0 0],[ax(3) ax(4)],'Color',linecolors(8,:),'LineWidth',3)
axis manual
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
        strvcat(printSnames{currcond_nums2==group}),...
        'Color',linecolors(group_ind,:),'FontWeight','Bold',...
        'FontSize',12)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'o','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'x','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
end
stddev2 = sqrt((std(scores(:,2))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
stddev1 = sqrt((std(scores(:,1))^2)*finv(out_cutoff,pcs,obs-pcs)*pcs*obs/(obs-pcs));
title([printDname ' ' typename ' Scores ('  ', ' num2str(totGenes)  ')'],'FontWeight','Bold',...
    'FontSize',24)
xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
x = [-stddev1:stddev1/50:stddev1];
yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
yl = -yu;
plot(x,yu,'Color',linecolors(8,:),'LineWidth',3)
plot(x,yl,'Color',linecolors(8,:),'LineWidth',3)
hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')

set(gca,'XTick',[],'YTick',[],'LineWidth',2)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs12_scores.pdf'])
currfig = [figinf '_pcs12_scores'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])

 
%% Loadings
% subplot(1,10,[7:10])


norms = sqrt(sum(coefs(:,1:2).^2,2));
[normsort sortinds] = sort(norms,'descend');
toPlot = sortinds(1:totWeights);

figure(2), clf
load_scale = min([ax(2)/max(coefs(coefs(:,1)>0,1))...
    ax(1)/min(coefs(coefs(:,1)<0,1))...
    ax(4)/max(coefs(coefs(:,2)>0,2)),...
    ax(3)/min(coefs(coefs(:,2)<0,2))]);     


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
text(load_scale*coefs(toPlot,1)',...
    load_scale*coefs(toPlot,2)',cellstr(strvcat(genesc{toPlot})),'FontWeight','Bold')
title([strrep(dataset,'_',' ') ' ' typename ' Loadings ('  ', ' num2str(totGenes)  ')'],'FontWeight','Bold',...
    'FontSize',24)

xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
currfig = [figinf '_pcs12_loadings'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])
 
%% Loadings 1, 2
% % subplot(1,10,[7:10])


norms = sqrt(sum(coefs(:,1:2).^2,2));
[normsort sortinds] = sort(norms,'descend');
toPlot = sortinds(1:totWeights);

figure(2), clf
load_scale = min([ax(2)/max(coefs(coefs(:,1)>0,1))...
    ax(1)/min(coefs(coefs(:,1)<0,1))...
    ax(4)/max(coefs(coefs(:,2)>0,2)),...
    ax(3)/min(coefs(coefs(:,2)<0,2))]);     


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
text(load_scale*coefs(toPlot,1)',...
    load_scale*coefs(toPlot,2)',cellstr(strvcat(genesc{toPlot})),'FontWeight','Bold')
title([strrep(dataset,'_',' ') ' ' typename ' Loadings ('  ', ' num2str(totGenes)  ')'],'FontWeight','Bold',...
    'FontSize',24)
xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
currfig = [figinf '_pcs12_loadings'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])


%% Scores Plot with Ellipse -- PCs 3 & 4
figure(3), clf
out_cutoff = .95;
% subplot(1,10,2:5)
%% PC plane 3, 4
% get axis set
axis auto
plot(1.1*min(scores(:,3)),1.2*min(scores(:,4)),'w')
hold on
plot(1.2*max(scores(:,3)),1.1*max(scores(:,4)),'w')        
ax = axis();
plot([ax(1) ax(2)],[0 0],'Color',linecolors(8,:),'LineWidth',3)
plot([0 0],[ax(3) ax(4)],'Color',linecolors(8,:),'LineWidth',3)
axis manual
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    text(scores(currcond_nums2==group,3),scores(currcond_nums2==group,4),...
        replace_underscore_w_spaces(cellstr(strvcat(printSnames{currcond_nums2==group}))),...
        'Color',linecolors(group_ind,:),'FontWeight','Bold',...
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
title([printDname ' ' typename ' Scores ('  ', ' num2str(totGenes)  ')'],'FontWeight','Bold',...
    'FontSize',24)
xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
x = [-stddev1:stddev1/50:stddev1];
yu = sqrt((stddev2^2)*(1-(x.^2)/(stddev1^2)));
yl = -yu;
plot(x,yu,'Color',linecolors(8,:),'LineWidth',3)
plot(x,yl,'Color',linecolors(8,:),'LineWidth',3)
hold on, plot(-[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot(+[stddev1 stddev1],[ax(3) ax(4)],'k:')
hold on, plot([ax(1) ax(2)],[stddev2 stddev2],'k:')
hold on, plot([ax(1) ax(2)],-[stddev2 stddev2],'k:')

set(gca,'XTick',[],'YTick',[],'LineWidth',2)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
currfig = [figinf '_pcs34_scores'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])

%% PC Plane 3/4 Loadings

totWeights = 60;

norms = sqrt(sum(coefs(:,3:4).^2,2));
[normsort sortinds] = sort(norms,'descend');
toPlot = sortinds(1:totWeights);

figure(3), clf
load_scale = min([ax(2)/max(coefs(coefs(:,3)>0,1))...
    ax(1)/min(coefs(coefs(:,3)<0,1))...
    ax(4)/max(coefs(coefs(:,4)>0,2)),...
    ax(3)/min(coefs(coefs(:,4)<0,2))]);     


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
text(load_scale*coefs(toPlot,3)',...
    load_scale*coefs(toPlot,4)',cellstr(strvcat(genesc{toPlot})),'FontWeight','Bold')
title([strrep(dataset,'_',' ') ' ' typename ' Loadings ('  ', ' num2str(totGenes)  ')'],'FontWeight','Bold',...
    'FontSize',24)
xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)
ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
    'FontSize',20)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [9 8]);
% print('-painters','-dpdf','-r1800',[figinf '_pcs34_scores.pdf'])
currfig = [figinf '_pcs34_loadings'];
save2pdf([currfig '.pdf'],gcf,150)
saveas(gcf,currfig,'fig')
disp([currfig '.pdf'])

%% PLOT 3 axis at once!

figure(200),clf
plot3(1.1*min(scores(:,1)),1.2*min(scores(:,2)),1.2*min(scores(:,3)),'w')
hold on
plot3(1.2*max(scores(:,1)),1.1*max(scores(:,2)),1.2*max(scores(:,3)),'w')    
plot3(scores(:,1),scores(:,2),scores(:,3),'w.')
for samp = 1: length(scores(:,1))
    plot3([0 scores(samp,1)],[0,scores(samp,2)],[0,scores(samp,3)],'y-')
end


% for ob = 1:obs
%     plot3([0 scores(ob,1)],[0 scores(ob,2)],[0 scores(ob,3)],'Color',[.75 .75 .75])
% end

hold on
for group_ind = 1:totgroups
    group = group_nums(group_ind);
    
    text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
        scores(currcond_nums2==group,3),...
        strvcat(printSnames{currcond_nums2==group}),...
        'Color',linecolors(group_ind,:),'FontWeight','Bold',...
        'FontSize',14)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'o','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
%     text(scores(currcond_nums2==group,1),scores(currcond_nums2==group,2),...
%         'x','FontWeight','Bold','Color',...
%         colors(group_ind,:),'FontSize',70)
end

grid on

ax = axis();
plot3([ax(1) ax(2)],[0 0],[0 0],'Color',linecolors(8,:),'LineWidth',3)
plot3([0 0],[ax(3) ax(4)],[0 0],'Color',linecolors(8,:),'LineWidth',3)
plot3([0 0],[0 0],[ax(5) ax(6)],'Color',linecolors(8,:),'LineWidth',3)

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
% title([printDname ' ' typename ' Loadings ('  ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
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
% % title([printDname ' ' typename ' Loadings ('  ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
% %     'FontSize',24)
% % xlabel(['PC 1 (' roundstring2(var_exp(1)) '%)'] ,'FontWeight','Bold',...
% %     'FontSize',20)
% % ylabel(['PC 2 (' roundstring2(var_exp(2)) '%)'] ,'FontWeight','Bold',...
% %     'FontSize',20)
% 
% in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');
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
% title([printDname ' ' typename ' Loadings ('  ', ' num2str(frac_measure_minimum)  ')'],'FontWeight','Bold',...
%     'FontSize',24)
% xlabel(['PC 3 (' roundstring2(var_exp(3)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% ylabel(['PC 4 (' roundstring2(var_exp(4)) '%)'] ,'FontWeight','Bold',...
%     'FontSize',20)
% 
% in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');
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
% %% Pareto Diagram -- Variance Explained
% % if plots(4)
% %     T = scores;
% %     P = coefs;
% %     denom = sum(diag(dataset*dataset'));
% %     var_exp = zeros(pcs,1);
% %     cum_var_exp = zeros(pcs,1);
% %     for pc = 1:pcs
% %         t = T(:,pc);
% %         p = P(:,pc);
% %         var_exp(pc) = sum(diag((t*p')*(t*p')'))/denom;
% %         if pc > 1
% %             cum_var_exp(pc) = var_exp(pc) + cum_var_exp(pc-1);
% %         else
% %             cum_var_exp(pc) = var_exp(pc);
% %         end
% %     end
% % %     subplot(3,6,15)
% %     figure(10)
% %     bar(100*var_exp)
% %     hold on
% %     plot(100*cum_var_exp,'bo')
% %     plot(100*cum_var_exp,'b')
% %     set(gca,'XTick',1:pcs,'XTickLabel',num2str([1:pcs]'),...
% %         'FontWeight','Bold')
% %     axis tight
% %     xlabel('Principle Components','FontWeight','Bold')
% %     ylabel('% Varience Explained','FontWeight','Bold')
% %     title('Variance Explained','FontWeight','Bold')
% % end
% 

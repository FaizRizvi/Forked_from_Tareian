%% visPriorOverlapsVsModelSize_BBLS_final_anyPrior

%% Debugging
clear all
% close all

%% Path to required functions
path2ems_functions = fullfile('/Users/caz3so/Dropbox/z_storage/rotations/Miraldi/scripts/emily_functions');

%% plotting parameters
regWidth = 3; % for p-r plots
controlWidth = 5;
fontSizeBar = 18; % for AUPR bar graphs
xSize = 8;       % for P-R plots
ySize = 4.75;
lineWidth = 3; 
figPause = 0;  % figPause --> 0, don't set figure for pdf, figPause --> 1, you can drag the figure windows to make the figure larger for pdf

darkBlue = [0 0 170]./255;
mediumBlue = [0 85 255]./255;
lightBlue = [0.5843    0.8157    0.9882];
darkRed = [170 0 0]./255;
mediumRed = [228 26 28]./255;
pink = [ 0.9686    0.5059    0.7490];
mediumGreen = [0 .45 0];
darkGreen = [0 .75 0];
lightGreen = [0 .25 0];

fontSize = 14;
%% Inf P-R results
maxModelSize = 25;
totGenes = 1999;

instCut = .1;  % cutoff for instability
outBit = [''];
timeOpt = ''; % added to get older versions of Th17 TRNs
addOutBit = '_PR';
kfoldCV = 5;
bsTot = 50;

manuscriptDir = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2';
totIntsInPrior = 5791;
disp('')

% params: (1) weight (g-prior or bias), (2) TFA option, (3) BBSR-BIC or
% LASSO-StARS, (4) line style, (5) line color, (6) line width
params = {1,'','bStARS','-',pink,regWidth;
    .5,'','bStARS','-',mediumRed,regWidth;
    .25,'','bStARS','-',darkRed,regWidth;
    .5,'_TFmRNA','bStARS','-',mediumBlue,regWidth;
    .25,'_TFmRNA','bStARS','-',lightBlue,regWidth;
    1,'_TFmRNA','bStARS','-',darkBlue,regWidth};

dataID = 'LO';
priorList = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/inputs/priors/priorlist.txt';
priorLocation = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/inputs/priors';
% priorList = ['/Users/emiraldi/erm/MariaP/Inferelator/input/priorLists/priorList_1709validationTEXT_CV.txt'];
%  = [1 5 6 7 8];% 11];%[1 2 5 8 11];%26];%:13 18:19];  % correspond to rows in the priorList
priorInds = [1];%:8];% 5];% 8];%[8];% 2];%[6:8]; % 8 --> ChIP, 1 --> sA(Th17) 1, 2, 6-8
% priorInds = [2:7 9 10 12];

tfCutoff = 10;
maxTfs = 20; % maximum TFs on x axis
%% Performance Analysis  
% add path to performance and other functions
addpath(path2ems_functions) % add path to functions
load linecolors.mat         % matrix of plotting colors for figures
load cbcolors
controlColor = .66*[1 1 1];

% /Users/emiraldi/erm/MariaP/Inferelator/input/lymph36broad_cv/Th17_TARGpap17v0FC0p58FDR10_REGSexp80nPapFC0p58FDR25_ivTh/Results_lassoStARS_ext_SS100/stabNetworksDebug_sG_15gene/ChBod10kMac10bin_bias10_sp.tsv
extraInf = '';
outInf = '';
resultFolder = '';
infModelDirLS = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/outputs/networks_targ0p05_SS50_GF_LI_LO_bS5';
noPriorResultsLS = '';
infModelDirBB = '';%'/Users/emiraldi/erm/MariaP/Inferelator/output/Th17_TARGpap17v0FC0p58FDR10_REGSexp80nPapFC0p58FDR25_ivTh';
noPriorResultsBB = '';%'Th17_48h_cut4_sA_p5_huA_TFmRNA_bs_50_weight_1_SS';
totalSSbb = 50;
outDir = fullfile(manuscriptDir,['ModelStats']);
mkdir(outDir)


%% read in prior file
pIn = fopen(priorList,'r');
C = textscan(pIn,'%s%s%s','Delimiter','\t');
fclose(pIn);
%%
priorFileTexts = C{1};
priorNames = C{2};
priorTsvs = C{3};
totPoi = length(priorFileTexts);
% priorFile = priorFileTexts{tind};
%%
for priorInd = priorInds
priorFile = priorFileTexts{priorInd};
priorName = priorNames{priorInd};   
priorTsv = priorTsvs{priorInd};
totInfResults = size(params,1);
% legendInf = cell(totInfResults,1);

%% open prior file and get interactions
% fid = fopen(fullfile(priorLocation,priorTsv),'r');
% % get first line and see what regulators we have    
% tline = fgetl(fid);    
% %pRegsTmp = strsplit(tline,'\t');
% % pRegsTmp = cellstr(strvcat(pRegsTmp)); % get rid of first \t if it exists
% totPRegs = length(strsplit(tline,'\t'));        
% fclose(fid);
% % get the rest of the data using textscan
% fid = fopen(fullfile(priorLocation,priorTsv),'r');
% C = textscan(fid,[repmat('%s',1,totPRegs)],'Delimiter','\t','Headerlines',1);
% fclose(fid);
% pRegs = C{1};
% pTargs = C{2};  
% priorInteractions = strcat(pRegs',pTargs');

    figure(1), clf % P-R
    hold on
    plot(tfCutoff*[1 1],[0 100],'--','Color',.83*[1 1 1],'LineWidth',3)
    legendInf = 'TF/gene_{cutoff}';
%     % plot random P-R ( : grey )
%     plot([0 1],gsInfs(gsInd).randPR*[1 1],':','LineWidth',lineWidth+1,'Color',controlColor)
%     hold on, grid on
%     set(gca,'FontSize',fontSize,'FontWeight','Bold')
%     xlabel('Recall','FontWeight','Bold','FontSize',fontSize)
%     ylabel('Precision','FontWeight','Bold','FontSize',fontSize)
%     axis([ 0 1 0 1])
%     grid on, grid minor
            
for tind = 1:totInfResults
    weight = params{tind,1};
    
    tfaOpt = params{tind,2};
    methodInf = params{tind,3};
    currMark = params{tind,4};
    currColor = params{tind,5};
    currWidth = params{tind,6};
    if and(find(ismember({'bStARS'},methodInf)),isnan(weight)) % LS no prior
        infModelDir = infModelDirLS;
        noPriorResults = noPriorResultsLS;
        infInput = fullfile(infModelDir,[noPriorResults extraInf '_sp.tsv']);
        disp([methodInf ''])
        fid = fopen(infInput,'r');
        C = textscan(fid,'%s%s%f','Delimiter','\t','Headerlines',1);
        regsNet = C{1};
        targsNet = C{2};
        quantilesRefinedTmp = C{3};
        totInts = length(quantilesRefinedTmp);
        currInts = strcat(regsNet',targsNet');
        inPriorVecTmp = ismember(currInts,priorInteractions);
        %% now need to find the prior... slow goin'                
        [vals, inds] = sort(abs(quantilesRefinedTmp),'descend');
        quantilesRefined = quantilesRefinedTmp(inds);
        inPriorVec = inPriorVecTmp(inds);

    elseif find(ismember({'bStARS'},methodInf)) % LS w/prior
        infModelDir = infModelDirLS;
        infInput = fullfile(infModelDir,[priorFile '_bias' num2str(100*weight) tfaOpt extraInf '.mat']);
        load(infInput)
        disp('LS')
%         gompers
%         fid = fopen(infInput,'r');
%         C = textscan(fid,'%s%s%f%f','Delimiter','\t','Headerlines',1);
%         quantilesRefinedTmp = C{3};
%         inPriorVecTmp = C{4};
%         [vals, inds] = sort(abs(quantilesRefinedTmp),'descend');
%         quantilesRefined = quantilesRefinedTmp(inds);
%         inPriorVec = inPriorVecTmp(inds);

    elseif and(find(ismember({'BB'},methodInf)),isnan(weight)) % BB no prior
        noPriorResults = fullfile(infModelDirBB,noPriorResultsBB);
        infModelDir = infModelDirBB;
        fullInfNet = ls(fullfile(noPriorResults,['*.tsv']));
        disp([methodInf ' No Prior'])
        fullInfNet = cellstr(strvcat(fullInfNet)); infInput = fullInfNet{1};  
        fid = fopen(infInput,'r');
        C = textscan(fid,[repmat('%s',1,19)],'Delimiter','\t','Headerlines',1);
        fclose(fid);    
%         disp('Fix this!')
%         inPriorVec = str2double(C{18});
        regsNet = C{1};
        targsNet = C{2};
        quantilesRefinedTmp = C{3};        
        currInts = strcat(regsNet',targsNet');
        inPriorVec = ismember(currInts,priorInteractions);

    elseif find(ismember({'BB'},methodInf)) % BB w/ prior
        fullInfNetGlob = fullfile(infModelDir,...
            [ priorFile tfaOpt '_bs_' num2str(totalSSbb) '_weight_' num2str(weight) '_SS'],['*.tsv']);
        disp(fullInfNetGlob)
        try
            fullInfNet = ls(fullInfNetGlob);
        catch
            error(['Check Inferelator model for ' priorName])
        end        
        fullInfNet = cellstr(strvcat(fullInfNet)); infInput = fullInfNet{1};
        fid = fopen(infInput,'r');
        C = textscan(fid,[repmat('%s',1,19)],'Delimiter','\t','Headerlines',1);
        fclose(fid);    
        inPriorVec = str2double(C{18});
    else
        error('Inference method not recognized as "BB" or "LS".')
    end

%     if find(ismember({'BB'},methodInf))
%         fid = fopen(infInput,'r');
%         C = textscan(fid,[repmat('%s',1,19)],'Delimiter','\t','Headerlines',1);
%         fclose(fid);    
%         inPriorVec = str2double(C{18});
% %         quantilesRefined = str2double(C{3});
%     else
%         fid = fopen(infInput,'r');
%         C = textscan(fid,'%s%s%f%f','Delimiter','\t','Headerlines',1);
%         quantilesRefinedTmp = C{3};
%         inPriorVecTmp = C{4};
%         [vals, inds] = sort(abs(quantilesRefinedTmp),'descend');
%         quantilesRefined = quantilesRefinedTmp(inds);
%         inPriorVec = inPriorVecTmp(inds);
%     end
    
    totEdges = length(inPriorVec);        
    currMaxMod = min(maxModelSize,floor(totEdges/totGenes)); % max TFs / gene model
    currModSizes = [1:currMaxMod]';
    percentPrior = zeros(currMaxMod,1);
    percentPriorRecall = zeros(currMaxMod,1);
%     percentNeg = zeros(currMaxMod,1);

    for ms = currModSizes'
        inds = 1:totGenes*ms;
        percentPrior(ms) = 100*length(find(inPriorVec(inds)))/(totGenes*ms);
        percentPriorRecall(ms) = 100*length(find(inPriorVec(inds)))/totIntsInPrior;
%         percentNeg(ms) = 100*length(find(quantilesRefined(inds)<0))/(totGenes*ms);
    end
        
    figure(1) % plot percent prior information
    plot(currModSizes,percentPrior,...
        currMark,'LineWidth',currWidth,'Color',currColor)
    plot(currModSizes,percentPriorRecall,...
        ':','LineWidth',currWidth,'Color',currColor)

    legendInf = strvcat(legendInf,['% Net in Prior' strrep(tfaOpt,'_',', ')],...
        ['% Prior in Net' strrep(tfaOpt,'_',', ')]);
    
end

    figure(1), % P-R
    title(strvcat([priorName strrep(extraInf,'_',' ')], ' '),'FontSize',fontSize+2)
    xlabel('Mean TFs / Gene','FontSize',fontSize+2,'FontWeight','Bold')
    ylabel('% Edges','FontSize',fontSize+2,'FontWeight','Bold')
    set(gca,'FontSize',fontSize,'FontWeight','Bold')
    plot(tfCutoff*[1 1],[200 201],'--','Color',.83*[1 1 1],'LineWidth',3)
    axis([0 maxTfs 0 100])
    grid on, grid minor, box on
    legendInf = strvcat(legendInf,'TF/gene_{cutoff}');
    legend(legendInf,'Location','EastOutside')
    currFile = fullfile(outDir,[priorFile '_priorPerModSize' extraInf]);
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
    print('-painters','-dpdf','-r150',[currFile '.pdf'])
    saveas(gcf,currFile,'fig')
    disp(currFile)
        
end
% end
disp('Finished')
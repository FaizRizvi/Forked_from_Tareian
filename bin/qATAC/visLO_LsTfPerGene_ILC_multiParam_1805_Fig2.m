%% visLO_LsTfPerGene_ILC_multiParam_1804
% visualize LO performance 
% Goal: compare models -- plot: 
%   1. R2^pred vs. network size 
%   2. SSE_pred vs. network size
% plot as subplots (unlabeled)

clear all
path_to_emsFXNs = fullfile('/Users/caz3so/Dropbox/z_storage/rotations/Miraldi/scripts/emily_functions');
addpath(path_to_emsFXNs)

% modelMode = 'BBSR-BIC'; % OPTIONS: 'BBSR-BIC' OR 'LASSO-StARS'
saveFig = 1;

kfoldCvOpt = [1];% 1] % 0 --> LO = exVivo Th, 1 --> CV = 5

% data location info
kfoldCV = 1;
bsTot = 50; % for BBSR-BIC
ssTot = 50; % subsamples for LASSO-StARS
ncuts = 50; % number of models bit based on bootstrap cutoffs
outInf = '';
maxTfs = 20;

%% plotting colors of code
load cbcolors.mat
darkBlue = [0 0 170]./255;
mediumBlue = [0 85 255]./255;
lightBlue = [0.5843    0.8157    0.9882];
darkRed = [170 0 0]./255;
mediumRed = [228 26 28]./255;
pink = [ 0.9686    0.5059    0.7490];
    % colors here go dark to light
blues = [darkBlue; mediumBlue; lightBlue];
reds = [darkRed; mediumRed; pink];
browns = [cbcolors(7,:)/2; cbcolors(7,:); mean([cbcolors(7,:);1 1 1])];
purples = [cbcolors(3,:)/2; cbcolors(3,:); mean([cbcolors(3,:);1 1 1])];
greens = [cbcolors(5,:)/2; cbcolors(5,:); mean([cbcolors(5,:);1 1 1])];
regWidth = 3; % for p-r plots
currWidth = 3;

manuscriptDir = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2';

instabText = '';

infModelDirBB = '';%'/Users/emiraldi/erm/MariaP/Inferelator/output/Th17_TARGpap17v0FC0p58FDR10_REGSexp80nPapFC0p58FDR25_ivTh';
noPriorResultsBB = '';%'Th17_48h_cut4_sA_p5_huA_TFmRNA/LOorCVname_CV1_bs_50_weight_1_SS_nCuts50.mat';
infModelDirLS = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/outputs/instabilities_targ0p05_SS50_GF_LI_LO_bS5';
totalSS = 50;
lambdaMax = 10;        % maximum lambda penalty starting place
lambdaMin = 1E-3;%1E-4;       % minimum lambda penalty starting place
meanEdgesPerGene = 20;
instabSource = 'Network';
targInstability = .05;
noPriorResultsLS = '';
priorList = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/inputs/priors/priorlist.txt';
infNetInds = [1];%:8];%26];%:13 18:19];  % correspond to rows in the priorList

% outBase = ['/Users/emiraldi/erm/MariaP/Inferelator/input/ILCs1707/ILCs_FDR10_FC1/infStARSfigs_CV' num2str(kfoldCV) dataID outInf];

% params: (1) weight (g-prior or bias), (2) TFA option, (3) BBSR-BIC or
% LASSO-StARS, (4) line style, (5) line color, (6) line width
params = {...    NaN,'','BB','-',lightBlue,currWidth; % NaN --> No Prior
    ...NaN,'','LS','-',pink,currWidth; % NaN --> No Prior
    1,'','bStARS','-',pink,regWidth;
    .5,'','bStARS','-',darkRed,regWidth;
    .25,'','bStARS','-',darkRed,regWidth;
    .5,'_TFmRNA','bStARS','-',mediumRed,regWidth;
    .25,'_TFmRNA','bStARS','-',mediumRed,regWidth;
    1,'_TFmRNA','bStARS','-',mediumBlue,regWidth};

tfCutoff = 10;

totParams = size(params,1);
paramInds = 1:totParams;

% plot parameters
fontSize = 14;
lineWidth = 2;
% xSize = 15;
% ySize = 2.5;
xSize = 8;       % for P-R plots
ySize = 4.75;

% ignore the second column
loFolders = {'GF_LI_LO','LI'};
            % 'GF_SI_LO','SI'
            % 'HH_LI_LO','HHLI'
            % 'SFB_SI_LO','SFBSI'};

%% end inputs
% mkdir(outBase)
%% read in prior file
pIn = fopen(priorList,'r');
C = textscan(pIn,'%s%s','Delimiter','\t');
fclose(pIn);
priorFiles = C{1};
priorNames = C{2};
totPoi = length(priorNames);

totInfResults = length(infNetInds);
legendInf = cell(totInfResults+1,1);
legendInf{1} = 'TF Cutoff';
% legendInf{2} = 'No Prior';
% need to get some info on priors to start and so will load first results

% totSubPlots = size(params,1);

markers = {'-',':','o-','-+'};

outDir = fullfile(manuscriptDir,['R2anal_Fig2']);
mkdir(outDir)
disp(outDir)

for nind = 1:totInfResults
    currInd = infNetInds(nind);
    priorFile = priorFiles{currInd};
    priorName = priorNames{currInd};
    currColor = cbcolors(rem(nind,length(cbcolors))+1,:);
    marker = markers{ceil(nind/length(cbcolors))};

    for loind = 1:size(loFolders,1)

        outCV = loFolders{loind,1};
        % limits = loFolders{loind,2};
        % minEdge = limits(1);
        % maxEdge = limits(2);
        % minSSE = limits(3);
        % maxSSE = limits(4);
        printLO = loFolders{loind,2};
        zoomAxisR2 = loFolders{loind,2};
        titleCV = outCV;    

        figure(1), clf, hold on % for R^2
        % figure(2), clf, hold on % for SSE
        
        plot(tfCutoff*[1 1],[0 1],'--','Color',.83*[1 1 1],'LineWidth',currWidth)

        for sp = paramInds
        weight = params{sp,1};
        tfaOpt = params{sp,2};
        infMeth = params{sp,3};
        lineStyle = params{sp,4};
        color = params{sp,5};
        currWidth = params{sp,6};

        % titleBase = [infMeth ', w=' num2str(weight) strrep(tfaOpt,'_',', ')];
        figBase = [priorFile '_' outCV ];

            toPlot = 1;
            if length(find(ismember({'bStARS'},infMeth))) > 0
                if isnan(weight) % No Prior
                    currMat = strrep(fullfile(infModelDirLS,noPriorResultsLS),'LOorCVname',outCV);
                    legendInf{sp+1} = ['No Prior'];
                else
                    currMat = fullfile(infModelDirLS,...
                        [priorFile '_bias' num2str(weight*100) ...
                            tfaOpt '_r2pred.mat']);
                    ls(currMat)
                    legendInf{sp+1} = [priorName ', bias = ' num2str(weight) strrep(tfaOpt,'_',', ')];

                end
            end
            load(currMat)
            totGenes = size(paramNumsMean,2);
            figure(1), 
            currR2pred = r2_pred;
            plot(sum(paramNumsMean,2)/totGenes,currR2pred,lineStyle,'Color',color,'LineWidth',currWidth)

        end

        % figure(2), % SSE_pred
        % % legend(legendInf,'Location','EastOutside')
        % axis([minEdge maxEdge minSSE maxSSE])
        figure(1), % R2_pred
%         plot(tfCutoff*[1 1],[20 21],'--','Color',.83*[1 1 1],'LineWidth',currWidth)
        hold on, grid on, set(gca,'FontSize',fontSize,'FontWeight','Bold'), grid minor, box on
        axis([0 maxTfs 0 1])
        xlabel('Mean TFs / Gene','FontSize',fontSize+2)
        ylabel('R^2_{pred}','FontSize',fontSize+2)
        title(strvcat([priorName ', LO = ' printLO], '   '),'FontSize',fontSize + 2)
        
%         legendInf = strvcat(strvcat(legendInf),'TF/gene','TF/gene');
        
        


%         plot(tfCutoff*[1 1],[20 21],'-','Color',.5*[1 1 1],'LineWidth',currWidth)
%         plot(tfCutoff*[1 1],[20 21],':','Color',.5*[1 1 1],'LineWidth',currWidth)
        
        legend(legendInf,'Location','EastOutside')
        if saveFig       
            figure(1)
            currFig = fullfile(outDir,[figBase '_r2pred']); 
            saveas(gcf,currFig,'fig')
            fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
            print('-painters','-dpdf','-r600',[currFig '.pdf'])
            disp(currFig)
            
            

%             figure(1)
%             for pp = paramInds
%                 subplot(length(paramInds),1,pp)
%                 axis(zoomAxisR2)
%             end
%             currFig = fullfile(outDir,[figBase '_r2pred_zoomAx']); 
%             saveas(gcf,currFig,'fig')
%             fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
%             print('-painters','-dpdf','-r150',[currFig '.pdf'])
%             disp(currFig)


        %     figure(2)
        %     currFig = fullfile(outDir,[figBase '_SSEpred']); 
        %     saveas(gcf,currFig,'fig')
        %     fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
        %     print('-painters','-dpdf','-r150',[currFig '.pdf'])
        %     disp(currFig)

        %     figure(2)
        %     
        %     currFig = fullfile(outDir,[figBase '_legend']); 
        %     saveas(gcf,currFig,'fig')
        %     fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
        %     print('-painters','-dpdf','-r150',[currFig '.pdf'])
        %     disp(currFig)
        end
        
    end
end
% figure(2)
% legend(legendInf,'Location','East')
% currFig = fullfile(outDir,[figBase '_LEGEND']); 
% saveas(gcf,currFig,'fig')
% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [10 20]);
% print('-painters','-dpdf','-r150',[currFig '.pdf'])
% disp(currFig)

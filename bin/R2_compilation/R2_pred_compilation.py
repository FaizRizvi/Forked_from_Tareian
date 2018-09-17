#!/usr/bin/env python

import pandas as pd 

"""" visualize LO performance """"
#   Goal: compare models -- plot: 
#   1. R2^pred vs. network size 
#   2. SSE_pred vs. network size
#   plot as subplots (unlabeled)

"""File Locations and Parameters"""
manuscriptDir = '/Users/caz3so/Dropbox/thesis/data/20180521_Myeloid_inference_sortingLO';
infModelDirLS = '/Users/caz3so/Dropbox/thesis/data/20180518_Myeloid_inference_R2/outputs/instabilities_targ0p05_SS50_SFB_SI_LO_bS5';
Prior_meta = '/Users/caz3so/Dropbox/thesis/data/20180521_Myeloid_inference_sortingLO/inputs/priors/Prior_meta.csv';
R2_meta = 'R2_meta.csv'
LO_meta = 'LO_meta.csv'

maxTfs = 20;
meanEdgesPerGene = 20;
tfCutoff = 10;

infNetInds = 1;  # correspond to rows in the priorList

# The parameters are defined by a metadata csv and are are: (1) weight (g-prior or bias), (2) TFA option, 
# (3) BBSR-BIC or LASSO-StARS, (4) line style, (5) line color, (6) line width
params = pd.read_csv(R2_meta, header = 0)

totParams = len(params,1);
paramInds = 1:totParams;

# plot parameters
fontSize = 14;
lineWidth = 2;
xSize = 8; 
ySize = 4.75;

loFolders = pd.read_csv(LO_meta, header=0)
           
"""Read in Prior File"""
pIn = pd.read_csv(Prior_meta, header=0)
priorFiles = pIn['priorfiles'];
priorNames = pIn['priornames'];
totPoi = len(priorNames);

totInfResults = len(infNetInds);
legendInf = cell(totInfResults+1,1);
legendInf{1} = 'TF Cutoff';

markers = {'-',':','o-','-+'};

outDir = fullfile(manuscriptDir,['R2anal_Fig2']);
mkdir(outDir)

for nind = 1:totInfResults
    currInd = infNetInds(nind);
    priorFile = priorFiles{currInd};
    priorName = priorNames{currInd};

    for loind = 1:size(loFolders,1)

        outCV = loFolders{loind,1};
        printLO = loFolders{loind,2};
        zoomAxisR2 = loFolders{loind,2};
        titleCV = outCV;    

        figure(1), clf, hold on 
        
        plot(tfCutoff*[1 1],[0 1],'--','Color',.83*[1 1 1],'LineWidth',currWidth)

        for sp = paramInds
        weight = params{sp,1};
        tfaOpt = params{sp,2};
        infMeth = params{sp,3};
        lineStyle = params{sp,4};
        color = params{sp,5};
        currWidth = params{sp,6};
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


        figure(1), % R2_pred
        
        hold on, grid on, set(gca,'FontSize',fontSize,'FontWeight','Bold'), grid minor, box on
        axis([0 maxTfs 0 1])
        xlabel('Mean TFs / Gene','FontSize',fontSize+2)
        ylabel('R^2_{pred}','FontSize',fontSize+2)
        title(strvcat([priorName ', LO = ' printLO], '   '),'FontSize',fontSize + 2)

        legend(legendInf,'Location','EastOutside')
        if saveFig       
            figure(1)
            currFig = fullfile(outDir,[figBase '_r2pred']); 
            saveas(gcf,currFig,'fig')
            fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
            print('-painters','-dpdf','-r600',[currFig '.pdf'])
            disp(currFig)
        end
        
    end
end

 
% loadGeneExprWithMetadata

clear all
close all

addpath(fullfile('~','Dropbox','Data','DC_GEO_Data','jake'))

%% Load quantitative data
inputGeneExpressionFile = 'KK221_expression.txt';
outputFolder = 'kk21';
dataset = 'kk21';  % add a dataset name, this will be used in figure titles and output files

% First see how many conditions / columns are in the file by reading the
% first line with tgetl
fid = fopen(inputGeneExpressionFile,'r');
tline = fgetl(fid);
fclose(fid);
conditionsc = strsplit(tline,'\t');    % get sample names
conditionsc = {conditionsc{2:end}}'; % leave out first field ("Name") as it doesn't refer to sample
totSamps = length(conditionsc); 
% Now that we know # of samples, reopen and get the rest
fid = fopen(inputGeneExpressionFile,'r');
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',1);
fclose(fid);
genesc = C{1};          % get gene names
ncounts = [C{2:end}];   % get gene expression values

%% Load metadata file and get cat names
listing = 'KK21_meta.txt';
fid = fopen(listing);
tline=fgetl(fid);
sampMetaCats = strsplit(tline,'\t');
totCats = length(sampMetaCats);
fclose(fid);

%% Open file and import data
fid = fopen(listing);
tmpSampMetaData = textscan(fid,['%s' repmat('%s',1,totCats)],'Delimiter','\t','Headerlines',1);
fclose(fid);
sampNames = tmpSampMetaData{1};
sampMetaData = cell(totSamps,totCats);
for samp = 1:totSamps
    currInd = find(ismember(conditionsc,sampNames{samp}));
    for cate = 1:totCats
        sampMetaData{currInd,cate} = tmpSampMetaData{cate}{samp};
    end
end

matName = fullfile(outputFolder,[dataset '_geneExprMeta']);

%% Save File as .mat
save([matName '.mat'],...
   'conditionsc',...                  11x1                   1544  cell                
   'dataset',...                       1x10                    20  char                
   'genesc',...
   'ncounts',...                  344941x11              30354808  double              
   'sampMetaCats',...                  1x8                   1044  cell                
   'sampMetaData')%,...                  1x9                   1132  cell    
disp([matName '.mat created!! :)'])

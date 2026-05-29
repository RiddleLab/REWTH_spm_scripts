function [outputBnormFiles] = backNormROI_REWTH(structuralImage,outputSpaceImage,toBackNormROIs,outputPrefix,outputDirectory,varargin)
% INPUTS:
%   structuralDirectory - directory for structural - generate norm files here
%   structuralImage - image that will be normalized
%   meanEpiFile - image for the mean of the functionals
%   preprocessingDirectory - make a directory for the backNormROIs
%   toBackNormROIs - images that will be backnormalized
%                   assume MNI spaces
%   outputPrefix - append this prefix to the backNormROIs
%
% OPTIONAL INPUTS:
%   restart = 1/0 default is 0, this will delete any existing ROIs
%
%   FLOW OF SCRIPTS      
%   1.  Normalize the structural image
%           generate norm and inverse norm parameters
%
%   2.  Back-normalize the Atlas ROIs into structural space
%
%   3.  Binarize the Atlas ROI in structural space
%
%   4.  Resample the ROI from Structural into EPI space
%

% OPTIONAL INPUTS
if length(varargin)>=1
    restart = varargin{1};
else
    restart = 0;
end

% Number of ROIs to be back-normalized
numROIs = length(toBackNormROIs);
outputBnormFiles = cell(1,numROIs);

% ROIs saved in sub directory within structural directory
[structuralDirectory,structFilename,~] = fileparts(structuralImage);
structuralDirectory = [structuralDirectory '/'];

%% Step 1: estimate the noramlization of the structural image
%   get the inverse normalization estimate as well
inverseDeformPrefix = 'iy_';
inverseDeformationFndr = dir([structuralDirectory inverseDeformPrefix structFilename '*']);
if isempty(inverseDeformationFndr)
    fprintf('Running Segmentation\n');
    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {structuralImage};
    % mean bias corrected image is written out
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    % change defaults to generate [inverse forward] normalization
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
    spm_jobman('run',matlabbatch);
    inverseDeformationFndr = dir([structuralDirectory inverseDeformPrefix structFilename '*']);
end
inverseDeformationFile = [structuralDirectory inverseDeformationFndr(1).name];

%% STEP 2: Apply inverse normalization
for roiIdx = 1:numROIs
    toBnormFile = toBackNormROIs{roiIdx};
    [~,fileName,ext]=fileparts(toBnormFile);
    outputFile = [outputDirectory outputPrefix fileName ext];
    outputBnormFiles{roiIdx} = outputFile;
    if restart && exist(outputFile,'file')==2
        delete(outputFile);
    end
    if exist(outputFile,'file')~=2
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {inverseDeformationFile}; %'/home/despo/riddler/PATRC/fMRI/Proc/sub07/Preproc_oct2016/Structural/iy_onssub07_patrc-0003-00001-000160-01.nii'};
        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {structuralImage}; %'/home/despo/riddler/PATRC/fMRI/Proc/sub07/Preproc_oct2016/Structural/monssub07_patrc-0003-00001-000160-01.nii'};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = {toBnormFile};
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {outputDirectory};%{'/home/despo/riddler/PATRC/fMRI/Proc/sub07/Preproc_oct2016/Structural/ROIs'};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {outputSpaceImage};%{'/home/despo/riddler/PATRC/fMRI/Proc/sub07/Preproc_oct2016/Structural/monssub07_patrc-0003-00001-000160-01.nii'};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        spm_jobman('run',matlabbatch);
        fprintf('Done\n')
        
        movefile([outputDirectory 'w' fileName ext],outputFile);
    end
end
end
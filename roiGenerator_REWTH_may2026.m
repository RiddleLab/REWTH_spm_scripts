%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate ROIs for TMS-targeting in REWTH
%   Justin Riddle & Timothy McDermott
%   Florida State University
%   Inquiries: justin.riddle@fsu.edu & tmcdermott@fsu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Goal:
%       Create regions of interest for rmPFC & dlPFC
%
%   Steps:
%       1. Manually type in the rmPFC coodinates in MNI space
%       2. Normalize structural scan and then back-normalize the rmPFC ROI
%       3. Calculate seed-based connectivity for sgACC
%       4. Manually type in the dlPFC coordinates in MNI space
%       5. Back-normalize the dlPFC site into subject space
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optional flags that change the data pipeline
% version control
% This will go on the folders and on the files themselves
glmVersion = 'april2026';
roiVersion = 'april2026';

%% Data paths
% Raw data is stored on CHAOS drive
RAW_DATA = 'R:\REWTH\';
% Scripts (where this script lives) should be in GitHub synced folder
GITHUB_CODEBASE = 'C:\Users\tjm25d\Documents\';
% Local data analysis on computer
LOCAL_DATA = 'C:\Users\tjm25d\Documents\REWTH_FMRI\REWTH_FMRI_SPM\';

% Load relevant toolboxes
RIDDLER_TOOLBOX = [GITHUB_CODEBASE 'RiddlerToolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));
SPM12_TOOLBOX = [GITHUB_CODEBASE 'Toolboxes/spm12/'];
addpath(SPM12_TOOLBOX);
spm('defaults', 'FMRI');
MARSBAR_TOOLBOX = [SPM12_TOOLBOX 'toolbox/marsbar/'];
addpath(MARSBAR_TOOLBOX);
marsbar('on');

% Raw data directories
BIDS_MRI = [RAW_DATA 'MRI_Data_BIDS\'];
RAW_BEHAV = [RAW_DATA 'Raw_Behavior\'];

% Output processed directories
GLM_MRI =[LOCAL_DATA 'GLM_MRI\'];
mkdir_JR(GLM_MRI);
% Directory for ROIs
ROI_DIR = [GLM_MRI 'ROIs_' glmVersion '/'];
mkdir_JR(ROI_DIR);
GENERAL_ROIS = [LOCAL_DATA 'GeneralROIs/'];
mkdir_JR(GENERAL_ROIS);

% Design matrices for the task
DES_MAT = [GLM_MRI 'DesignMatrices\'];
mkdir_JR(DES_MAT);

% Subjects
SUBJECTS = {'sub006','sub002', 'sub007','sub013'};
numSub = length(SUBJECTS);

% Two different types of analysis
TYPES_NAMES = {'Savor','Rest'};
numTypes = length(TYPES_NAMES);

validAnswer = 0;
while ~validAnswer
    inputSubID_str = input('Enter the participant subID number (e.g., 7): ','s');
    inputSubID = str2double(inputSubID_str);
    if ~isnan(inputSubID) && (inputSubID > 0 && inputSubID < 1000)
        validAnswer = 1;
    else
        validAnswer = 0;
    end
end

% Subject code
subID = sprintf('%03d',inputSubID);
subject = ['sub' subID];

validAnswer = 0;
while ~validAnswer
    inputREDO_str = input(sprintf('Do you want to redo ROIs for %s: ',subject),'s');
    if strcmpi(inputREDO_str,'y')
        validAnswer = 1;
        redoROIs_FLAG = 1;
    elseif strcmpi(inputREDO_str,'n')
        validAnswer = 1;
        redoROIs_FLAG = 0;
    else
        validAnswer = 0;
    end
end

% Directory to save out the ROIs
SUB_ROI = [ROI_DIR subject '/'];
mkdir_JR(SUB_ROI);

% BIDS fMRI prep - paths to preprocessed data
SUB_PREPROC_MRI = [BIDS_MRI 'BIDS_' subject '\derivatives\fmriprep-25.2.0\sub-' subID '\'];
SUB_PREPROC_FUNC = [SUB_PREPROC_MRI 'func\'];
SUB_PREPROC_STRUCT = [SUB_PREPROC_MRI 'anat\'];

mniSpaceStr = 'MNI152NLin2009cAsym_res-2';

% Participant's MNI anatomical
subMNI_structFilename = sprintf('sub-%s_space-%s_desc-preproc_T1w.nii',...
    subID,mniSpaceStr);
subMNI_structFile = [SUB_ROI subMNI_structFilename];
chaos_structFile = [SUB_PREPROC_STRUCT subMNI_structFilename];
gzip_subMNI_structFile = [SUB_PREPROC_STRUCT subMNI_structFilename '.gz'];
if exist(subMNI_structFile,'file')~=2
    if exist(chaos_structFile,'file')~=2 && exist(gzip_subMNI_structFile,'file')==2
        gunzip(gzip_subMNI_structFile);
    end
    movefile(chaos_structFile,subMNI_structFile);
end

% Participant's subject-space anatomical
subSpace_structFilename = sprintf('sub-%s_desc-preproc_T1w.nii',subID);
subSpace_structFile = [SUB_ROI subSpace_structFilename];
chaos_structFile = [SUB_PREPROC_STRUCT subSpace_structFilename];
gzip_subSpace_structFile = [SUB_PREPROC_STRUCT subSpace_structFilename '.gz'];
if exist(subSpace_structFile,'file')~=2
    if exist(chaos_structFile,'file')~=2 && exist(gzip_subSpace_structFile,'file')==2
        gunzip(gzip_subSpace_structFile);
    end
    movefile(chaos_structFile,subSpace_structFile);
end


%% Step 1. Manually type in the rmPFC coodinates in MNI space
% Then run beta series connectivity analysis
roiName = 'rmPFC';
coords = {'X','Y','Z'};
bNorm_rmPFC_roiFile = [SUB_ROI subject '_subjectSpace_' roiName '_' roiVersion '.nii'];
mni_rmPFC_roiFile = [SUB_ROI subject '_mniSpace_' roiName '_' roiVersion '.nii'];
if redoROIs_FLAG && exist(bNorm_rmPFC_roiFile,'file')==2
    delete(bNorm_rmPFC_roiFile);
end
if redoROIs_FLAG && exist(mni_rmPFC_roiFile,'file')==2
    delete(mni_rmPFC_roiFile);
end
if exist(bNorm_rmPFC_roiFile,'file')~=2

    validAnswer = 0;
    while ~validAnswer
        inputROI_str = input(sprintf('Do you have %s coordinates for %s? (y or n): ',...
            roiName,subject),'s');
        if strcmpi(inputROI_str,'y')
            run_rmPFCcreation = 1;
            validAnswer = 1;
        elseif strcmpi(inputROI_str,'n')
            run_rmPFCcreation = 0;
            validAnswer = 1;
        end
    end

    % If the user is ready with the aMFG coordinates
    if run_rmPFCcreation
        % If it does not exist then generate it
        roiCoordinates = NaN(1,3);
        % Loop through each coordinate
        for coordIdx = 1:length(coords)
            coord = coords{coordIdx};
            % Request the user to type in the coordinates
            validAnswer = 0;
            while ~validAnswer
                thisCoord_str = input(sprintf('What is the %s for %s?: ',...
                    coord,roiName),'s');
                thisCoord = str2double(thisCoord_str);
                if ~isnan(thisCoord) && (abs(thisCoord) < 200)
                    roiCoordinates(coordIdx) = thisCoord;
                    validAnswer = 1;
                end
            end
        end
        roi_rootFilename = [subject '_MNI_' roiVersion];
        [outputImageFiles,~] = makeROIs_fromCoordinates(...
            roi_rootFilename,{'rmPFC'},roiCoordinates,'sphere',5,subMNI_structFile,SUB_ROI);
        defaultName_rmPFC_roiFile = outputImageFiles{1};
        movefile(defaultName_rmPFC_roiFile, mni_rmPFC_roiFile)

        %% Step 2. Normalize structural scan and then back-normalize the rmPFC ROI
        % files for structural get written into the STRUCT_DIR 'ROIs/'
        % epis go into reslice dir
        [backNormROIfiles] = backNormROI_REWTH(subSpace_structFile,...
            subSpace_structFile,{mni_rmPFC_roiFile},'bNorm_',SUB_ROI);
        movefile(backNormROIfiles{1},bNorm_rmPFC_roiFile);

    end
end

%% Step 3: check that the MNI subgenual cingulate exists
sgACC_MNI_roiFile = [GENERAL_ROIS 'mniSpace_sgACC_may2026.nii'];
if exist(sgACC_MNI_roiFile,'file')~=2
    % In Cash et al., 2021 & in Fox et al., 2012 the peak
    % coordinate in MNI space was (6, 16, -10) with a 10 mm sphere
    roiName = 'sgACC';
    % If it does not exist then generate it
    roiCoordinates = NaN(1,3);
    % Loop through each coordinate
    for coordIdx = 1:length(coords)
        coord = coords{coordIdx};
        % Request the user to type in the coordinates
        validAnswer = 0;
        while ~validAnswer
            thisCoord_str = input(sprintf('What is the %s for %s?: ',...
                coord,roiName),'s');
            thisCoord = str2double(thisCoord_str);
            if ~isnan(thisCoord) && (abs(thisCoord) < 200)
                roiCoordinates(coordIdx) = thisCoord;
                validAnswer = 1;
            end
        end
    end
    [outputImageFiles,~] = makeROIs_fromCoordinates(...
        'MNI',{'sgACC'},roiCoordinates,'sphere',10,subMNI_structFile,GENERAL_ROIS);
    movefile(outputImageFiles{1},sgACC_MNI_roiFile);
end

%% Step 4: Segment the MNI T1w before creating a nuissance mask
numTissues = 5;
allTissuesExist = 1;
for tissueIdx = 1:numTissues
    if exist(sprintf('%sc%i%s',SUB_ROI,tissueIdx,subMNI_structFilename),'file')~=2
        allTissuesExist = 0;
    end
end
if ~allTissuesExist
    segmentNormalize_spm12(subMNI_structFile);
end

%% Step 5: Create a nuissance mask for this participant in MNI space
mni_sub_brainMaskFile = [SUB_ROI 'mni_nuissanceMask_' subject '_' roiVersion '.nii'];
if exist(mni_sub_brainMaskFile,'file')~=2
    segmentedTissueFiles = cell(1,numTissues);
    for tissueIdx = 1:numTissues
        segmentedTissueFiles{tissueIdx} = ...
            sprintf('%sc%i%s',SUB_ROI,tissueIdx,subMNI_structFilename);
    end
    % GM, WM, CSF, SKULL, OUTSIDE
    weights = [0.2 -0.9 -0.9 -0.9 -0.9];
    % Generate a nuissance mask
    nuissanceMaskMaker(weights,segmentedTissueFiles,mni_sub_brainMaskFile);
end

%% Step 6: Run Seed based connectivity with sgACC
SUB_GLM = [GLM_MRI 'Rest_MNI_' glmVersion '/' subject '/'];
mni_restingStateResiduals = [SUB_GLM 'Residuals_4d.nii'];
sgACC_seedConnFile = [SUB_ROI 'mni_sgACC_seedConn_Rest_' subject '_' roiVersion '.nii'];
if exist(sgACC_seedConnFile,'file')~=2
    connectivitySeedBased({sgACC_MNI_roiFile}, {mni_restingStateResiduals},mni_sub_brainMaskFile,{sgACC_seedConnFile})
end

%% Step 7: Manually type in the dlPFC coordinates in MNI space
roiName = 'dlPFC';
coords = {'X','Y','Z'};
bNorm_dlPFC_roiFile = [SUB_ROI subject '_subjectSpace_' roiName '_' roiVersion '.nii'];
if redoROIs_FLAG && exist(bNorm_dlPFC_roiFile,'file')==2
    delete(bNorm_dlPFC_roiFile);
end

mni_dlPFC_roiFile = [SUB_ROI subject '_mniSpace_' roiName '_' roiVersion '.nii'];
if redoROIs_FLAG && exist(mni_dlPFC_roiFile,'file')==2
    delete(mni_dlPFC_roiFile);
end
if exist(bNorm_dlPFC_roiFile,'file')~=2
    validAnswer = 0;
    while ~validAnswer
        inputROI_str = input(sprintf('Do you have %s MNI coordinates for %s? (y or n): ',...
            roiName,subject),'s');
        if strcmpi(inputROI_str,'y')
            run_dlPFCcreation = 1;
            validAnswer = 1;
        elseif strcmpi(inputROI_str,'n')
            run_dlPFCcreation = 0;
            validAnswer = 1;
        end
    end
end

% If the user is ready with the aMFG coordinates
if run_dlPFCcreation
    % If it does not exist then generate it
    roiCoordinates = NaN(1,3);
    % Loop through each coordinate
    for coordIdx = 1:length(coords)
        coord = coords{coordIdx};
        % Request the user to type in the coordinates
        validAnswer = 0;
        while ~validAnswer
            thisCoord_str = input(sprintf('What is the %s MNI coord for %s?: ',...
                coord,roiName),'s');
            thisCoord = str2double(thisCoord_str);
            if ~isnan(thisCoord) && (abs(thisCoord) < 200)
                roiCoordinates(coordIdx) = thisCoord;
                validAnswer = 1;
            end
        end
    end
    roi_rootFilename = [subject '_MNI_' roiVersion];
    [outputImageFiles,~] = makeROIs_fromCoordinates(...
        roi_rootFilename,{'dlPFC'},roiCoordinates,'sphere',5,subMNI_structFile,SUB_ROI);
    defaultName_dlPFC_roiFile = outputImageFiles{1};
    movefile(defaultName_dlPFC_roiFile,mni_dlPFC_roiFile)

    %% Step 8: back-normalize the dlPFC ROI
    % files for structural get written into the STRUCT_DIR 'ROIs/'
    % epis go into reslice dir
    [backNormROIfiles] = backNormROI_REWTH(subSpace_structFile,...
        subSpace_structFile,{mni_dlPFC_roiFile},'bNorm_',SUB_ROI);
    movefile(backNormROIfiles{1},bNorm_dlPFC_roiFile);

    %% Step 9: run seed-based connectivity from the dlPFC ROI
    % And visually confirm that the connectivity pattern matches
    % canonical frontal-parietal control network
    dlPFC_seedConnFile = [SUB_ROI 'mni_dlPFC_seedConn_Rest_' subject '_' roiVersion '.nii'];
    if redoROIs_FLAG && exist(dlPFC_seedConnFile,'file')==2
        delete(dlPFC_seedConnFile);
    end
    if exist(dlPFC_seedConnFile,'file')~=2
        connectivitySeedBased({mni_dlPFC_roiFile}, {mni_restingStateResiduals},mni_sub_brainMaskFile,{dlPFC_seedConnFile})
    end
end % if running dlPFC ROI creation (coords acquired)
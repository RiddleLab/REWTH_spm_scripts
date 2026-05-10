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
RAW_DATA = 'Z:\REWTH\';
% Scripts (where this script lives) should be in GitHub synced folder
GITHUB_CODEBASE = 'C:\Users\jr23z\GitHub_Codebase\';
% Local data analysis on computer
LOCAL_DATA = 'C:\Users\jr23z\LocalDataAnalysis\REWTH\';

% Load relevant toolboxes
RIDDLER_TOOLBOX = [GITHUB_CODEBASE 'RiddlerToolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));
SPM12_TOOLBOX = [GITHUB_CODEBASE 'Toolboxes/spm12/'];
addpath(SPM12_TOOLBOX);
spm('defaults', 'FMRI');
MARSBAR_TOOLBOX = [GITHUB_CODEBASE 'Toolboxes/marsbar-0.44/'];
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

% Design matrices for the task
DES_MAT = [GLM_MRI 'DesignMatrices\'];
mkdir_JR(DES_MAT);

% Subjects
SUBJECTS = {'sub006'};
numSub = length(SUBJECTS);

% Two different types of analysis
TYPES_NAMES = {'Savor','Rest'};
numTypes = length(TYPES_NAMES);

% Loop through each subject
for subIdx = 1:numSub
    subject = SUBJECTS{subIdx};
    subID = subject(4:end);
    
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

    % Participant's MNI-space functional

    %% Step 1. Manually type in the rmPFC coodinates in MNI space
    % Then run beta series connectivity analysis
    roiName = 'rmPFC';
    coords = {'X','Y','Z'};

    bNorm_rmPFC_roiFile = [SUB_ROI subject '_subjectSpace_' roiName '_' roiVersion '.nii'];
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
            rmPFC_roiFile = outputImageFiles{1};
        
            %% Step 2. Normalize structural scan and then back-normalize the rmPFC ROI
            % files for structural get written into the STRUCT_DIR 'ROIs/' 
            % epis go into reslice dir
            [backNormROIfiles] = backNormROI_REWTH(subSpace_structFile,...
                subSpace_structFile,{rmPFC_roiFile},'bNorm_',SUB_ROI);
            movefile(backNormROIfiles{1},bNorm_rmPFC_roiFile);
            
        end

        % 3. Calculate seed-based connectivity for sgACC
        % 4. Manually type in the dlPFC coordinates in MNI space
        % 5. Back-normalize the dlPFC site into subject space
    end
end
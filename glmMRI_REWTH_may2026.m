%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLM for REWTH MRI
%   Justin Riddle & Timothy McDermott
%   Florida State University
%   Inquiries: justin.riddle@fsu.edu & tmcdermott@fsu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Goal:
%       Preprocessing for fMRI has been run through fMRIPREP
%       Now run the general linear model and task contrasts
%
%   Next steps:
%       Plot contrast and select ROIs
%           rmPFC from savor vs view: hedonic
%           Draw ROI in subgenual cingulate
%           Then from sgACC calculate peak anti-correlation in aMFG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optional flags that change the data pipeline
% version control
% This will go on the folders and on the files themselves
glmVersion = 'april2026';

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

% Raw data directories
BIDS_MRI = [RAW_DATA 'MRI_Data_BIDS\'];
RAW_BEHAV = [RAW_DATA 'Raw_Behavior\'];

% Output processed directories
GLM_MRI =[LOCAL_DATA 'GLM_MRI\'];
mkdir_JR(GLM_MRI);

% Design matrices for the task
DES_MAT = [GLM_MRI 'DesignMatrices\'];
mkdir_JR(DES_MAT);

% Subjects
SUBJECTS = {'sub002', 'sub006', 'sub007','sub013'};
numSub = length(SUBJECTS);

% Two different types of analysis
TYPES_NAMES = {'Savor','Rest'};
numTypes = length(TYPES_NAMES);
TYPES_TASKNAME = {'savor','rest_run-0'};
TYPES_NUMRUNS = [3 2];

% Task parameters
SAVOR_COND = {'Savor','View'};
numSavor = length(SAVOR_COND);
VALENCE_COND = {'Eudaimonic','Hedonic','Neutral'};
numValence = length(VALENCE_COND);

% Assemble all conditions
ALL_COND = cell(1,numSavor*numValence);
condIdx = 0;
for savorIdx = 1:numSavor
    savor = SAVOR_COND{savorIdx};
    for valenceIdx = 1:numValence
        valence = VALENCE_COND{valenceIdx};
        condIdx = condIdx + 1;
        ALL_COND{condIdx} = [savor valence];
    end
end
numCond = length(ALL_COND);

% framewise displacement threshold - flag these
FrameDisp_threshold = 0.2;

% Loop through each subject
for subIdx = 1:numSub
    subject = SUBJECTS{subIdx};
    subID = subject(4:end);
    if strcmp(subject, 'pilot002') %%%% temporary to process Jake pilot data - TM on 5-18-26
        subID = 'pilot002';
    else
        subID = subject(4:end);
    end

    % BIDS fMRI prep - paths to preprocessed data
    SUB_PREPROC_MRI = [BIDS_MRI 'BIDS_' subject '\derivatives\fmriprep-25.2.0\sub-' subID '\'];
    SUB_PREPROC_FUNC = [SUB_PREPROC_MRI 'func\'];
    SUB_PREPROC_STRUCT = [SUB_PREPROC_MRI 'anat\'];

    % Subject raw behavioral files
    SUB_RAW_BEHAV = [RAW_BEHAV 'REWTH_' subject '\'];
    assert(exist(SUB_RAW_BEHAV,'dir')==7);

    for typeIdx = 1:numTypes
        runType = TYPES_NAMES{typeIdx};
        file_taskName = TYPES_TASKNAME{typeIdx};
        numRunsMax = TYPES_NUMRUNS(typeIdx);

        fprintf('\n\nRunning GLM for %s %s\n',subject,runType);

        % Output directories for this script
        GLM_DIR = [GLM_MRI runType '_MNI_' glmVersion '/'];
        mkdir_JR(GLM_DIR);

        % Subject specific output directories
        SUB_GLMDIR = [GLM_DIR subject '/'];
        mkdir_JR(SUB_GLMDIR);

        % The exhaustive analysis contains all of the condition of interest
        % Note that there is more statistical power for the neutral cue
        SUB_DESMAT = [DES_MAT  subject '\'];
        mkdir_JR(SUB_DESMAT);

        mriSpaceStr = 'MNI152NLin2009cAsym_res-2';
        % subject space version: mriSpaceStr = 'T1w';

        % Loop through the maximum number of EPIs and gather .nii files
        epiFiles = cell(1,0);
        numEPIs = 0;
        for epiIdx = 1:numRunsMax
            epiFile = sprintf('%ssub-%s_task-%s%i_space-%s_desc-preproc_bold.nii',...
                SUB_PREPROC_FUNC,subID,file_taskName,epiIdx,mriSpaceStr);
            gzip_epiFile = [epiFile '.gz'];
            if exist(epiFile,'file')~=2 && exist(gzip_epiFile,'file')==2
                fprintf('Unzipping MRI files for run %i\n\n',epiIdx);
                gunzip(gzip_epiFile);
            end
            if exist(epiFile,'file')==2
                numEPIs = numEPIs + 1;
                epiFiles{epiIdx} = epiFile;
            else
                % Break the for-loop if missing a run
                break; 
                % If there are weird recording issues, like runs are
                % miscounted, then this step will not work as intended
            end
        end

        %% Estimate the design matrix for the savor task or specify rest
        if strcmp(runType,'Savor')
            % Then make a design matrix
            % otherwise, the rest does not need this step
        
            % Design matrix files
            designMatrixFiles = cell(1,numEPIs);

            % Loop through each task run
            rawBehavFiles = cell(1,numEPIs);
            for epiIdx = 1:numEPIs
                % Raw behavior files
                rawBehavFndr = dir(sprintf('%sresults_SAVOR_%s_MRI_block%i_*.mat',...
                    SUB_RAW_BEHAV,subject,epiIdx));
                rawBehavFile = [SUB_RAW_BEHAV rawBehavFndr(1).name];
                rawBehavFiles{epiIdx} = rawBehavFile;

                % Output design matrix file
                designMatrixFile = sprintf('%s%s_designMatrix_savor_run%i_%s.mat',...
                    SUB_DESMAT,subject,epiIdx,glmVersion);
                designMatrixFiles{epiIdx} = designMatrixFile;

                % Make the design matrix if it does not exist
                if exist(designMatrixFile,'file')~=2

                    % Populate the onsets and durations for each condition
                    names = ALL_COND;
                    onsets = cell(1,numCond);
                    durations = cell(1,numCond);

                    % Load the results file
                    behavStruct = load(rawBehavFile);

                    % The time that the experiment started in computer time
                    experimentStartTIME = behavStruct.exptStart;

                    % Six "trials" that we call "mini-blocks"
                    numTrials = 6;
                    for trialIdx = 1:numTrials
    
                        % Extract relevant time points
                        trialInfo = behavStruct.(sprintf('miniBlock%i',trialIdx));
                        savorCond = trialInfo.savor;
                        valenceCond = trialInfo.valence;

                        condIdx = find(strcmpi(ALL_COND,[savorCond valenceCond]),1);

                        % The onset of this event is the onset of the first
                        % stimulus that is present
                        firstStimOnsetTIME = trialInfo.stimTiming(1,1);
                        % the offset of this event is the offeset of the
                        % last stimulus that is presneted (number 7)
                        lastStimOffsetTIME = trialInfo.stimTiming(end,end);

                        % Duration is from the first stim on to last stim off
                        trialDuration = lastStimOffsetTIME - firstStimOnsetTIME;
                        % Trial onset is calculated relative to expt start
                        trialOnset = firstStimOnsetTIME - experimentStartTIME;

                        % If it is the first instance
                        if isempty(onsets{condIdx})
                            % then make it this onset
                            onsets{condIdx} = trialOnset;
                            durations{condIdx} = trialDuration;
                        else
                            % otherwise append to the end
                            onsets{condIdx} = [onsets{condIdx} trialOnset];
                            durations{condIdx} = [durations{condIdx} trialDuration];
                        end
                    end
                    % save the design matrix file
                    desMatStruct = struct(...
                        'names',{names},...
                        'onsets',{onsets},...
                        'durations',{durations});
                    if exist(designMatrixFile,'file')==2
                        delete(designMatrixFile);
                    end
                    save(designMatrixFile,'-struct','desMatStruct');
                end % if the design matrix files does not exist
            end % loop through runs
           
            % Run the GLM based on these task maps
            extraContrasts = {...
                'AllCond_SavorEudaimonic_1_SavorHedonic_1_SavorNeutral_1_ViewEudaimonic_1_ViewHedonic_1_ViewNeutral_1',...
                'SavorEffect_SavorEudaimonic_1_SavorHedonic_1_SavorNeutral_1_ViewEudaimonic_-1_ViewHedonic_-1_ViewNeutral_-1',...
                'SavorEffectHedonic_SavorHedonic_1_ViewHedonic_-1',...
                'SavorEffectHedControlled_SavorHedonic_1_SavorNeutral_-1_ViewHedonic_-1_ViewNeutral_1',...
                'Extremes_SavorHedonic_1_ViewNeutral_-1',...
                'EudaimonicEffect_SavorEudaimonic_1_ViewEudaimonic_-1',...
                'HedonicEffect_SavorHedonic_1_SavorNeutral_-1_ViewHedonic_1_ViewNeutral_-1'};
            glmInfo = struct(...
                'extraContrasts',{extraContrasts},...
                'designMatrices',{designMatrixFiles},...
                'restingState',0);
        else
            glmInfo = struct(...
                'restingState',1);
        end

        %% Generate Nuissance regressors files
        regressorFiles = cell(1,numEPIs);
        numFD_rejection = zeros(1,numEPIs);
        % TM added line below in modified version on 5-10-26 to account for
        % different number of non-steady state outliers (3, 2, 1, or 0 - in order)
        allNumNuissanceRegressors = zeros(1,numEPIs);
        fdRejectionFile = [SUB_DESMAT subject '_' runType '_MNIspace_fdRejection_' glmVersion '.mat'];
        for epiIdx = 1:numEPIs
            % Write out a regressors file for this run
            glm_regressorsFile = sprintf('%sregressors_%s_%s_run%i_%s.txt',...
                SUB_DESMAT,subject,runType,epiIdx,glmVersion);
            if exist(glm_regressorsFile,'file')==2
                delete(glm_regressorsFile);
            end
            regressorFiles{epiIdx} = glm_regressorsFile;

            % Load up the nuissance regressors from fMRI prep
            fmriprep_regressorsFile = sprintf(...
                '%ssub-%s_task-%s%i_desc-confounds_timeseries.tsv',...
                SUB_PREPROC_FUNC,subID,file_taskName,epiIdx);
            fmriprep_regressorsFID = fopen(fmriprep_regressorsFile,'r');
            hdrLine = fgetl(fmriprep_regressorsFID);
            hdrParts = regexp(hdrLine,'\t','split');

            NUISSANCE_REGRESSORS = {...
                'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z',...
                'csf','white_matter','framewise_displacement'};

            % Sometimes there are non-steady state outliers that should be
            % added as additional regresssors
            nonSteady_hdrIdxs = find(contains(hdrParts,'non_steady_state_outlier'));
            numNonSteady = length(nonSteady_hdrIdxs);
            for nonSteadyIdx = 1:numNonSteady
                nonSteady_hdrIdx = nonSteady_hdrIdxs(nonSteadyIdx);
                % extract the header name
                regressorName = hdrParts{nonSteady_hdrIdx};
                % append to the end
                NUISSANCE_REGRESSORS{end+1} = regressorName;
            end

            % Let's try without global signal, the task data looks a lot better 'global_signal',
            numNuissanceRegressors = length(NUISSANCE_REGRESSORS);
            allNumNuissanceRegressors(epiIdx) = numNuissanceRegressors; % TM added on 5-10-26

            % First figure out which header is each regressor
            regressor_hdrIdxs = NaN(1,numNuissanceRegressors);
            for regIdx = 1:numNuissanceRegressors
                regressor = NUISSANCE_REGRESSORS{regIdx};
                regressor_hdrIdx = find(strcmp(regressor,hdrParts));
                regressor_hdrIdxs(regIdx) = regressor_hdrIdx;
            end

            % Then loop through and create the GLM regressor files with a
            % subset of the columns
            glm_regressorsFID = fopen(glm_regressorsFile,'w');
            moreLines = 1;
            firstLineFLAG = 1;
            while moreLines
                nextLine = fgetl(fmriprep_regressorsFID);
                if nextLine == -1
                    moreLines = 0;
                else
                    if firstLineFLAG
                        firstLineFLAG = 0;
                    else
                        fprintf(glm_regressorsFID,'\n');
                    end
                    fmriprep_parts = regexp(nextLine,'\t','split');
                    for regIdx = 1:numNuissanceRegressors
                        regressor = NUISSANCE_REGRESSORS{regIdx};
                        thisRegressor = str2double(fmriprep_parts{regressor_hdrIdxs(regIdx)});
                        if strcmp(regressor,'framewise_displacement')
                            thisRegressor = thisRegressor > FrameDisp_threshold;
                            if thisRegressor == 1
                                numFD_rejection(epiIdx) = numFD_rejection(epiIdx) + 1;
                            end
                        end

                        fprintf(glm_regressorsFID,'%0.20f',thisRegressor);
                        if (regIdx ~= numNuissanceRegressors)
                            fprintf(glm_regressorsFID,'\t');
                        end
                    end
                end
            end
            % Delete the final "enter" that goes to a blank line
            fclose(glm_regressorsFID);
            fclose(fmriprep_regressorsFID);
        end

        % Save out how many volumes are rejected from framewise
        % displacement
        if exist(fdRejectionFile,'file')==2
            delete(fdRejectionFile);
        end
        save(fdRejectionFile,'numFD_rejection');

        glmInfo.regressorFiles = regressorFiles;
        glmInfo.TR = 2;
        glmInfo.numSlices = 62; % for 62 slices with multiband 2
        glmInfo.refSlice = 1;
        if strcmp(runType,'Savor')
            glmInfo.numVolumes = 211;
        else
            glmInfo.numVolumes = 140;
        end
        glmInfo.numRegressors = allNumNuissanceRegressors;

        %% Run the GLM

        % If no REST scans, then skip this GLM and move on to SAVOR - % TM added on 5-18-26
        if numEPIs == 0
            fprintf('No EPI files found for %s %s - skipping to next GLM.\n', subject, runType);
            continue; % skip to next runType
        end

        glm_SPM12_JR_TM(SUB_GLMDIR, epiFiles, glmInfo);
    end % loop run types (task or rest)
end % subjects loop
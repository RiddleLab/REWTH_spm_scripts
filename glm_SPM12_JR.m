function [varargout] = glm_SPM12_JR(RESULTS_DIR, epiFiles, glmInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CREDITS
%
%       Justin Riddle
%       Florida State University
%       Last updated on April 20, 2026
%       Email for contact:  justin.riddle@fsu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RUN GENERAL LINEAR MODEL FOR FMRI
%
%   1.  SPECIFY FIRST LEVEL & MODEL ESTIMATION
%
%   2.  CONTRAST GENERATOR - one T-contrast for each condition in design (default)
%           (OPTIONAL create extra contrasts)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%preproc_SPM12_JR Run full preprocessing pipeline
%   Inputs:
%       outputDirectory - directory to create preprocessed data
%       glmInfo ~ includes the following options
%           &&&&&&& RESULTS PROCESSING &&&&&&&
%                  .designMatrices - {{1xnumEPIs}}
%                       paths to design matrices
%                       must include this to perform first level statistics
%                       onsets specified using seconds from first volume
%                           to view the finalized design matrix: SPM.xX.X
%                 ----ALTERNATIVE to glmInfo.designMatrices----
%                  .restingState = 1/0
%                       this will run resting state analysis - gets residuals
%                       no conditions here, just regressors
%                       turns on: segmentation regressors
%                                 includes "rest" in glmInfo.resultsName
%                                 changes spm12 directory to:
%                                   /home/despo/riddler/Toolboxes/spm_toolbox/spm12/
%                                   which keeps residuals instead of the
%                                   default to delete them
%                  .extraContrasts = cell string {1xnumContrasts}
%                       string format: contrastName_conditionName_value_conditionName_value
%                       example t-contrast: scenesMinusFaces_faces_-1_scenes_1
%                       example f-contrast: scenesByFaces_faces_-1_+_scenes_1
%                  .restartContrasts = 1/0  redo contrasts
%                       default = 0
%                  .regressorsFiles = cell string of paths to extra regressors text files
%                                           must have regressors equal to number of scans
%                  .restartResults = 1/0
%                       if 1, restart results but not image preprocessing
%                       default = 0
%                  .resultsName - string
%                       provide an alternative name for your results for
%                       multiple GLMs on the same data
%                       default = ''
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of EPI scans
numEPIs = length(epiFiles);

% Information required for timing
TR = glmInfo.TR;
numSlices = glmInfo.numSlices;
refSlice = glmInfo.refSlice;
numVolumes = glmInfo.numVolumes;
numRegressors = glmInfo.numRegressors;
if isscalar(numVolumes)
    allNumVolumes = ones(1,numEPIs)*numVolumes;
else
    allNumVolumes = numVolumes;
end

% Design matrices for task data (otherwise should be rest)
if isfield(glmInfo,'designMatrices')
    ALL_DES_MAT = glmInfo.designMatrices;
    restingStateFlag = 0;
else
    ALL_DES_MAT = {};
    restingStateFlag = 1;
end
% If resting state flag is 1, then run GLM without design matrices
% store the residuals for use afterwards
if isfield(glmInfo,'restingState')
    assert(restingStateFlag == glmInfo.restingState);
end
% This generates extra contrasts based on main effects
if isfield(glmInfo,'extraContrasts')
    extraContrasts = glmInfo.extraContrasts;
else
    extraContrasts = {};
end

% deletes old contrasts
if isfield(glmInfo,'restartContrasts')
    restartContrasts = glmInfo.restartContrasts;
else
    restartContrasts = 0;
end

% If you want to add extra regressors of your own do this with this variable
assert(isfield(glmInfo,'regressorFiles'))
regressorFiles = glmInfo.regressorFiles;

% Change mask threshold
if isfield(glmInfo,'resultsMaskThreshold')
    resultsMaskThreshold = glmInfo.resultsMaskThreshold;
else
    % this is the SPM default (very conservative)
    resultsMaskThreshold = 0.8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify First Level %%
%%%%%%%%%%%%%%%%%%%%%%%%%
glmFl = [RESULTS_DIR 'SPM.mat'];
if exist(glmFl,'file')
    fprintf('Skipping Specify First Level\n');
else
    fprintf('Running Specify First Level\n');
    clear matlabbatch
    spm_jobman('initcfg');
    for epiIdx = 1:numEPIs
        epiFile = epiFiles{epiIdx};
        [EPI_DIR, epiFilename, epiExt] = fileparts(epiFile);
        functionalName = [epiFilename epiExt];
        numVolumes = allNumVolumes(epiIdx);
        % collect final version of preprocessed epi files
        toSpecify = cell(numVolumes,1);
        for volIdx = 1:numVolumes
            toSpecify{volIdx} = [EPI_DIR '\' functionalName ',' num2str(volIdx)];
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(epiIdx).scans = toSpecify;

        % there might not be any design matrices for resting state scans
        if ~isempty(ALL_DES_MAT)
            matlabbatch{1}.spm.stats.fmri_spec.sess(epiIdx).multi = ALL_DES_MAT(epiIdx);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(epiIdx).multi_reg = regressorFiles(epiIdx);
    end
    % details of this specify first level call
    matlabbatch{1}.spm.stats.fmri_spec.dir = {RESULTS_DIR};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = numSlices;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = refSlice;
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = resultsMaskThreshold;
    spm_jobman('run',matlabbatch);
end

% must be at least one beta estimate for each regressor and EPI
if isempty(dir([RESULTS_DIR 'beta*']))
    %%%%%%%%%%%%%%%%%%%%%%
    %% MODEL ESTIMATION %%
    %%%%%%%%%%%%%%%%%%%%%%
    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {glmFl};
    if restingStateFlag
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1;
    end
    spm_jobman('run',matlabbatch);
    % this step will generate residual files for resting state
end

% if resting state
if restingStateFlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONVERT RESIDUALS TIMESERIES TO 4D %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    residuals4dFile = fullfile(RESULTS_DIR, 'Residuals_4d.nii');

    if ~isfile(residuals4dFile)
        fprintf('Converting residuals files into 4D timeseries...\n');

        % Find all residual images
        residualFileFndr = dir(fullfile(RESULTS_DIR, 'Res_*.nii'));

        % Full paths as cell array
        residualFiles = fullfile({residualFileFndr.folder}, {residualFileFndr.name});

        % Merge 3D residuals into one 4D NIfTI
        spm_file_merge(char(residualFiles), residuals4dFile);

        % Delete original 3D residual files
        delete(fullfile(RESULTS_DIR, 'Res_*.nii'));
    end
else
    % delete residuals if not running resting state
    residualPattern = [RESULTS_DIR 'ResI*.nii'];
    if ~isempty(dir(residualPattern))
        delete(residualPattern);
    end
end

% the following 3 steps all require conditions of interest
% not applicable for resting state without any events
% events for resting state could be a regressor - check for extra con
if (~isempty(ALL_DES_MAT) || ~isempty(extraContrasts))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BETA IMAGES - calculate name of each %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make a simple T-contrast for each condition versus baseline
    desMatNames = cell(1,0);
    % First gather all the unique names of conditions
    uniqueDesMatCond = 0;

    % this variable has all names of beta images in order
    namesInOrder = cell(1,1);
    orderIdx = 0;
    % loop through epis - betas in order by epi
    for epiIdx = 1:numEPIs
        if ~isempty(ALL_DES_MAT)
            % load design matrix for this epi
            DES_MAT = load(ALL_DES_MAT{epiIdx});
            % gather the condition names for novelty search
            thisDesMatNames = DES_MAT.names;
            thisNumCond = length(thisDesMatNames);
            % search through names in design matrix
            for desMatNameIdx = 1:thisNumCond
                thisDesMatName = thisDesMatNames{desMatNameIdx};
                orderIdx = orderIdx + 1;
                namesInOrder{orderIdx} = thisDesMatName;
                % gather novel conditions for contrast making
                if isempty(find(strcmpi(thisDesMatName,desMatNames),1))
                    % found a novel condition
                    uniqueDesMatCond = uniqueDesMatCond + 1;
                    % update novel names found in design matrices
                    desMatNames{uniqueDesMatCond} = thisDesMatName;
                end
            end
        end
        % each regressor gets a beta for each epi
        for regressIdx = 1:numRegressors
            orderIdx = orderIdx + 1;
            namesInOrder{orderIdx} = ['regressor' num2str(regressIdx)];
        end
    end
    % the last beta images are regressors for each epi
    for epiIdx = 1:numEPIs
        orderIdx = orderIdx + 1;
        namesInOrder{orderIdx} = ['regressEPI' num2str(epiIdx)];
    end
    % use these for contrast generation
    mainEffectContrasts = desMatNames;
    betaNamesInOrder = namesInOrder;
    numBetas = length(namesInOrder);

    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Contrast Generator %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numExtraContrasts = length(extraContrasts);
    numCond = length(mainEffectContrasts);
    numContrasts = numCond + numExtraContrasts;
    % get the names of every contrast
    contrastNames = mainEffectContrasts;
    for extraIdx = 1:numExtraContrasts
        extraContrastSpecs = extraContrasts{extraIdx};
        extraContrastSpecs = regexp(extraContrastSpecs,'_','split');
        contrastNames{numCond + extraIdx} = extraContrastSpecs{1};
    end
    % check that there is a contrast image for every main effect + extra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numExistingContrasts = length(dir([RESULTS_DIR 'con_*.nii']));
    if (numExistingContrasts >= numContrasts && ~restartContrasts)
        fprintf('Skipping Contrast Generator\n');
    else
        fprintf('Running Contrast Generator\n');
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.con.spmmat = {glmFl};
        for condIdx = 1:numCond
            condName = mainEffectContrasts{condIdx};
            % zero for each regressor between EPIs, zeroes at end equal to number of EPIs
            masterContrast = double(strcmpi(condName,betaNamesInOrder));
            matlabbatch{1}.spm.stats.con.consess{condIdx}.tcon.name = condName;
            matlabbatch{1}.spm.stats.con.consess{condIdx}.tcon.weights = masterContrast;
            matlabbatch{1}.spm.stats.con.consess{condIdx}.tcon.sessrep = 'none';
        end
        % Loop through extra contrasts
        for extraIdx = 1:numExtraContrasts
            masterContrast = zeros(1,numBetas);
            extraContrastSpecs = extraContrasts{extraIdx};
            extraContrastSpecs = regexp(extraContrastSpecs,'_','split');
            aspectIdx = 1; % skip name

            % The first value is the name of the contrast
            contrastName = extraContrastSpecs{aspectIdx};

            % Followed by a condition name, then a numeric value
            % and repeat for as many conditions as are in the contrast
            moreAspects = 1;
            while moreAspects
                aspectIdx = aspectIdx + 1;
                condName = extraContrastSpecs{aspectIdx};
                aspectIdx = aspectIdx + 1;
                aspectValue = str2double(extraContrastSpecs{aspectIdx});
                masterContrast = masterContrast ...
                    + (double(strcmpi(condName,betaNamesInOrder)) * aspectValue);

                if (aspectIdx+1) > length(extraContrastSpecs)
                    moreAspects = 0;
                end
            end
            matlabbatch{1}.spm.stats.con.consess{numCond + extraIdx}.tcon.name = contrastName;
            matlabbatch{1}.spm.stats.con.consess{numCond + extraIdx}.tcon.weights = masterContrast;
            matlabbatch{1}.spm.stats.con.consess{numCond + extraIdx}.tcon.sessrep = 'none';
        end
        % Delete the current contrasts
        matlabbatch{1}.spm.stats.con.delete = 1;
        spm_jobman('run',matlabbatch);
    end
end
end % end of function
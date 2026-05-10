%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second Level GLM for REWTH MRI - Savor Task
%   Timothy McDermott (with help from Claude)
%   Florida State University
%   Inquiries: tmcdermott@fsu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Goal:
%       Run second level (group) one-sample t-tests for each contrast
%       from the first level Savor task GLM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optional flags that change the data pipeline
% version control
glmVersion = 'april2026';

%% Subjects
SUBJECTS = {'sub006', 'sub002', 'sub007'};
numSub = length(SUBJECTS);

%% Data paths
GITHUB_CODEBASE = 'C:\Users\tjm25d\Documents\';
LOCAL_DATA = 'C:\Users\tjm25d\Documents\REWTH_FMRI\REWTH_FMRI_SPM\';

% Load relevant toolboxes
RIDDLER_TOOLBOX = [GITHUB_CODEBASE 'RiddlerToolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));
SPM12_TOOLBOX = [GITHUB_CODEBASE 'Toolboxes/spm12/'];
addpath(SPM12_TOOLBOX);
spm('defaults', 'FMRI');

% First level results directory
GLM_MRI = [LOCAL_DATA 'GLM_MRI\'];
FIRST_LEVEL_DIR = [GLM_MRI 'Savor_MNI_' glmVersion '\'];

% Second level output directory
SECOND_LEVEL_DIR = [GLM_MRI 'SecondLevel_Savor_MNI_' glmVersion '_N' num2str(numSub) '\'];
mkdir_JR(SECOND_LEVEL_DIR);

%% Contrast names - must match order in first level SPM.mat
CONTRAST_NAMES = {...
    'SavorEudaimonic',...
    'SavorHedonic',...
    'SavorNeutral',...
    'ViewEudaimonic',...
    'ViewHedonic',...
    'ViewNeutral',...
    'AllCond',...
    'SavorEffect',...
    'SavorEffectHedonic',...
    'SavorEffectHedControlled',...
    'Extremes',...
    'EudaimonicEffect',...
    'HedonicEffect'};
numContrasts = length(CONTRAST_NAMES);

%% Loop through each contrast
for conIdx = 1:numContrasts
    contrastName = CONTRAST_NAMES{conIdx};

    fprintf('\n\nRunning second level t-test for %s (N=%i)\n', contrastName, numSub);

    % Output directory for this contrast
    CON_DIR = [SECOND_LEVEL_DIR contrastName '\'];
    mkdir_JR(CON_DIR);

    % Skip if SPM.mat already exists
    if exist([CON_DIR 'SPM.mat'], 'file')
        fprintf('Skipping %s - already exists\n', contrastName);
        continue;
    end

    % Gather contrast images from each subject
    conImages = cell(numSub, 1);
    for subIdx = 1:numSub
        subject = SUBJECTS{subIdx};
        conImages{subIdx} = sprintf('%s%s\\con_%04d.nii', ...
            FIRST_LEVEL_DIR, subject, conIdx);
        assert(exist(conImages{subIdx}, 'file') == 2, ...
            'Missing contrast image for %s: %s', subject, conImages{subIdx});
    end

    %% Specify second level model (one-sample t-test)
    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.stats.factorial_design.dir = {CON_DIR};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conImages;
    spm_jobman('run', matlabbatch);

    %% Model estimation
    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[CON_DIR 'SPM.mat']};
    spm_jobman('run', matlabbatch);

    %% Contrast - mean activation across subjects
    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.stats.con.spmmat = {[CON_DIR 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [contrastName '_N' num2str(numSub)];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);

end % contrast loop
addpath('/projects/AnalysisTools/BCT/2014_04_05 BCT/')
subjNums = '013 014 016 017 018 021 023 024 026 027 028 030 031 032 033 034 035 037 038 039 040 042 043 045 046 047 048 049 050 053 055 056 057 058 062 063 066 067 068 069 070 072 074 075 076 077 078 079 081 085 086 087 088 090 092 093 094 095 097 098 099 101 102 103 104 105 106 108 109 110 111 112 114 117 119 120 121 122 123 124 125 126 127 128 130 131 134 135 137 138 139';
%subjNums = '013';
subjNumStr = strread(subjNums, '%s', 'delimiter', ' ');

for i=1:91
    subj=subjNumStr{i};
    subj
% This function runs a GLM on the Gordon Parcels (333) on a single subject
% This will only regress out noise parameters and HRF convolved unique 64 task set 
%
% Original purpose: To run activity-based MVPA on mean activity on ROIs within CCN Networks for each rule dimension.

% Input parameter:
%   subj = subject number as a string
%   gsr = 0 if no GSR, 1 if you want to include GSR 
%   nproc = number of processes to run in parallel -- will run (i.e., regress) either regions or voxels/vertices in parallel
    basedir = ['/projects/IndivRITL/data/'];
    analysisdir = [basedir '/' subj '/analysis/'];
    datadir = [basedir subj '/MNINonLinear/Results'];
    numTasks = 8;
    numTRs = numTasks*581;
    numROIs = 264;
    disp('Loading in and downsampling volume data to Power ROIs')
    disp('Loading Petersen lab 264 ROI set');
    %Using Freesurfer's NIFTI file handling
    psetHDR3D=load_nifti(['/projects/AnalysisTools/atlases/PetersenLab264_MNI_WithLabels.nii.gz']);
    %Reshape to 1D
    psetvol=reshape(psetHDR3D.vol,[size(psetHDR3D.vol,1)*size(psetHDR3D.vol,2)*size(psetHDR3D.vol,3) 1]);
    
    disp(['Loading task preprocessed fMRI data for subject ' subj]);
        dset4D=load_nifti(['/projects/IndivRITL/data/' subj '/analysis/Task_allruns.nii.gz']);
        %Reshape to 2D (space x time)
        dsetvol=reshape(dset4D.vol,[size(dset4D.vol,1)*size(dset4D.vol,2)*size(dset4D.vol,3) size(dset4D.vol,4)]);
        clear dset4D

        %Get ROI data
        ROIData_corr=zeros(numROIs,size(dsetvol,2));
        for regionNum=1:numROIs
            %Get region's data
            ROIData_corr(regionNum,:)=mean(dsetvol(psetvol==regionNum,:),1);
        end
        clear dsetvol
    
    disp(['Loading rest preprocessed fMRI data for subject ' subj]);
        dset4D=load_nifti(['/projects/IndivRITL/data/' subj '/MNINonLinear/Results/Rest1/Rest1.nii.gz']);
        %Reshape to 2D (space x time)
        dsetvol=reshape(dset4D.vol,[size(dset4D.vol,1)*size(dset4D.vol,2)*size(dset4D.vol,3) size(dset4D.vol,4)]);
        clear dset4D

        %Get ROI data
        ROIData_corr2=zeros(numROIs,size(dsetvol,2));
        for regionNum=1:numROIs
            %Get region's data
            ROIData_corr2(regionNum,:)=mean(dsetvol(psetvol==regionNum,:),1);
        end
        clear dsetvol
        
    dtseries_task = ROIData_corr;
    dtseries_rest = ROIData_corr2;
    
    %% Rest & Task GLM

    % Load only noise regressors and task regressors for data
    X = loadStimFiles_by64TaskSet_dhs264(subj,gsr);

    %import censor file
    TaskCensorfile=(['/projects/IndivRITL/data/motionfiles/' subj '_Task_censor.1D']);
    RestCensorfile=(['/projects/IndivRITL/data/motionfiles/' subj '_Rest_censor.1D']);
    TaskTRstokeep = importdata(TaskCensorfile);
    RestTRstokeep = importdata(RestCensorfile);
    TaskTRstokeep = logical(TaskTRstokeep);
    RestTRstokeep = logical(RestTRstokeep);
    
    %for loop to trim regressors
    for regNum=1:80
        newReg=X.regressors(:,regNum);
        trimReg=newReg(TaskTRstokeep==1);
        X.regressors_trim(:,regNum)=trimReg;
    end
    
    for restregNum=1:16
        newrestreg=X.noiseRegressors(:,restregNum);
        trimrestReg=newrestreg(RestTRstokeep==1);
        X.noiseRegressors_trim(:,restregNum)=trimrestReg;
    end
        
    %Now motion censoring on the timeseries
    for roiNum=1:264
        dtseries_task_reg(1,:)=dtseries_task(roiNum,:);
        dtseries_task_cen(roiNum,:)=dtseries_task_reg(TaskTRstokeep'==1);
    end
   
    for roiNum=1:264
        dtseries_rest_reg(1,:)=dtseries_rest(roiNum,:);
        dtseries_rest_cen(roiNum,:)=dtseries_rest_reg(RestTRstokeep'==1);
    end
    
    disp(['Running regression for subject ' subj]);
    
    % Instantiate empty arrays
    rest_resids = zeros(size(dtseries_rest_cen));
    rest_betas = zeros(numROIs, size(X.noiseRegressors_trim,2)+1);

    task_resids = zeros(size(dtseries_task_cen));
    disp(['Number of regressors in the task data is: ' num2str(size(X.regressors_trim,2))])
    task_betas = zeros(numROIs, size(X.regressors_trim,2)+1); 
    % Need to add one for constant regressor

    disp(['Begin regression for ' subj])
    
    % Begin for loop
    for regionNumb=1:264
        % Get the region's data
        ROITimeseries_task = dtseries_task_cen(regionNumb,:);
        ROITimeseries_rest = dtseries_rest_cen(regionNumb,:);
        % Regress out the task, keep the residuals, betas, and rsquare
        stats_task = regstats(ROITimeseries_task', X.regressors_trim, 'linear', {'r', 'beta', 'rsquare'});
        
        % Regress out noise parameters for resting state data
        stats_rest = regstats(ROITimeseries_rest', X.noiseRegressors_trim, 'linear', {'r', 'beta', 'rsquare'});

 
        % Collect rest regression results
        rest_resids(regionNumb,:) = stats_rest.r';
        rest_betas(regionNumb,:) = stats_rest.beta;

        % Collect task regression results
        task_resids(regionNumb,:) = stats_task.r';
        task_betas(regionNumb,:) = stats_task.beta';
    end
    
    disp(['FC calculation for ' subj])
    
    interactionMatrix_corrs_task=zeros(264,264,64);
    interactionMatrix_corrs_rest=zeros(264,264,1);
    interactionMatrix_pvals_task=zeros(264,264,64);
    interactionMatrix_pvals_rest=zeros(264,264,1);
    
    Trim_task_regs=X.regressors_trim(:,17:80);
    regressorMatrix_canonicalHRF_ByTask_binary=Trim_task_regs>0;
    
    for condNum=1:64
        taskTimePoints=task_resids(:,regressorMatrix_canonicalHRF_ByTask_binary(:,condNum));
        [interactionMatrix_corrs_task(:,:,condNum), interactionMatrix_pvals_task(:,:,condNum)] = corrcoef(taskTimePoints');
    end
    
    [interactionMatrix_corrs_rest(:,:,1), interactionMatrix_pvals_rest(:,:,1)] = corrcoef(rest_resids');
    
    disp(['Output results for ' subj])
    
    output.rest_betas=rest_betas;
    output.rest_resids=rest_resids;
    output.restFC_corr=interactionMatrix_corrs_rest;
    output.interactionMatrix_pvals_rest=interactionMatrix_pvals_rest;
    output.task_betas=task_betas;
    output.task_resids=task_resids;
    output.taskFC_corr=interactionMatrix_corrs_task;
    output.interactionMatrix_pvals_task=interactionMatrix_pvals_task;
    output.X=X;
    output.TaskTRstokeep=TaskTRstokeep;
    output.RestTRstokeep=RestTRstokeep;
    
    out_struct{i} = output;
    
    clear dtseries_rest
    clear dtseries_task
    clear dtseries_rest_cen
    clear dtseries_rest_reg
    clear dtseries_task_cen
    clear dtseries_task_reg
    clear X
    clear rest_resids
    clear task_resids
    clear regressorMatrix_canonicalHRF_ByTask_binary
    clear ROITimeseries_rest
    clear ROITimeseries_task
    clear taskTimePoints
    clear Trim_task_regs
    clear trimReg
    clear trimrestReg
    
end




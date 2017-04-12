function output = GordonROIGLM_64taskset_IndivRITL(subj, gsr, nproc)
% This function runs a GLM on the Gordon Parcels (333) on a single subject
% This will only regress out noise parameters and HRF convolved unique 64 task set 
%
% Original purpose: To run activity-based MVPA on mean activity on ROIs within CCN Networks for each rule dimension.

% Input parameter:
%   subj = subject number as a string
%   gsr = 0 if no GSR, 1 if you want to include GSR 
%   nproc = number of processes to run in parallel -- will run (i.e., regress) either regions or voxels/vertices in parallel
    
    numTasks = 8;
    numTRs = numTasks*581;

    disp('Loading in and downsampling surface data to Gordon Parcels')
    alldata = loadSurfaceData2(subj);
    dtseries_task = alldata.dtseries_task;
    dtseries_rest = alldata.dtseries_rest;
    numROIs = size(dtseries_task,1);
    
    %% Rest & Task GLM

    % Load only noise regressors and task regressors for data
    X = loadStimFiles_by64TaskSet_mc16(subj,gsr);

    % Instantiate empty arrays
    rest_resids = zeros(size(dtseries_rest));
    rest_betas = zeros(numROIs, size(X.noiseRegressors,2)+1);

    task_resids = zeros(size(dtseries_task));
    disp(['Number of regressors in the task data is: ' num2str(size(X.regressors,2))])
    task_betas = zeros(numROIs, size(X.regressors,2)+1); % Need to add one for constant regressor

    % Begin for loop
    parfor (regionNum=1:numROIs,nproc)
        % Get the region's data
        ROITimeseries_task = dtseries_task(regionNum,:);
        ROITimeseries_rest = dtseries_rest(regionNum,:);
        % Regress out the task, keep the residuals, betas, and rsquare
        stats_task = regstats(ROITimeseries_task', X.regressors, 'linear', {'r', 'beta', 'rsquare'});
        
        % Regress out noise parameters for resting state data
        stats_rest = regstats(ROITimeseries_rest', X.noiseRegressors, 'linear', {'r', 'beta', 'rsquare'});

 
        % Collect rest regression results
        rest_resids(regionNum, :) = stats_rest.r';
        rest_betas(regionNum,:) = stats_rest.beta;

        % Collect task regression results
        task_resids(regionNum,:) = stats_task.r';
        task_betas(regionNum,:) = stats_task.beta';
    end

    output.rest_resids = rest_resids;
    output.rest_betas = rest_betas;
    output.task_resids = task_resids;
    output.task_betas = task_betas;
end




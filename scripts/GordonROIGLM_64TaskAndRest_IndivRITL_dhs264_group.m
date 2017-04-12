%% Taku Ito
%% 3/13/15
%
%% Group analysis for modalityControl
%% Run GLM on all subjects and do some summary statistics
%
addpath('/projects/AnalysisTools/BCT/2014_04_05 BCT/')
subjNums = '013 014 016 017 018 021 023 024 025 026 027 028 030 031 032 033 034 035 037 038 039 040 042 043 045 046 047 048 049 050 052 053 055 056 057 058 062 062 064 066 067 068 069 070 072 074 075 076 077 078 079 081 085 086 087 088 090 092 093 094 095 097 098 099 101 102 103 104 105 106 108 109 110 111 112 114 115 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 134 135 136 137 138 139 140 141';
subjNumStr = strread(subjNums, '%s', 'delimiter', ' ');

out = cell(length(subjNumStr),1);

step1 = 1;
step2 = 1;
step3 = 0;
step4 = 0;
plot = 0;
writeSubjectDataToCSV=0;
allVerticesGLM_subjectwise=0;

numROIs=264;

% Runs subject-wise GLM for noise regression (rest) and task-evoked regression (task) in parallel
if step1 == 1
    %% Run in parallel
    for (i=1:length(out)) %Use maximum 5 processes
        out_nogsr{i} = GordonROIGLM_64TaskAndRest_IndivRITL_dhs264(subjNumStr{i},gsr, dtype, hipitchOnly, randomStims, taskRegression, treatRestWithTaskRegressors, blockdemeaning, excludeTaskEvokedResponseRegs);
    end
end

% Extracts the subject-wise correlation maps for both rest and task data (for each entire timeseries), and gets an average correlation map
if step2 == 1
    % Extract subject-wise data and get an average timeseries
    %corr_gsr = zeros(numROIs,numROIs,length(subjNumStr));
    corr_nogsr = zeros(numROIs,numROIs,length(subjNumStr));
    rest_nogsr = zeros(numROIs,numROIs,length(subjNumStr));

    for subj=1:length(subjNumStr)
        if covariance == 1
            corr_nogsr(:,:,subj) = cov(out_nogsr{subj}.regionData_taskRegressed_byMod');
            rest_nogsr(:,:,subj) = cov(out_nogsr{subj}.regionData_restRegressed');
        else
            %corr_nogsr(:,:,subj) = corrcoef(out_nogsr{subj}.regionData_taskRegressed_byMod');
            rest_nogsr(:,:,subj) = corrcoef(out_nogsr{subj}.regionData_restRegressed');
        end
        
    end

    %avg_corrMat = mean(corr_nogsr,3);
    avg_restMat = mean(rest_nogsr,3);

end

% Imports the parcellation network association, and prepares data for plotting (and re-orders parcels by network association)
if step3 == 1
    if strmatch(dtype,'Surface')
        % import parcels and reorder according to networks!!
        parcelfile = '/projects/ModalityControl/data/indivritl/Parcels/Parcels.xlsx';
        parcelimport = importdata(parcelfile);
        parcelimport = parcelimport.textdata;
        parcellabels = parcelimport(2:end,5);
        
        if strmatch(roi_order,'Power')
            out = reorderGordonToPower();
            xticks = out.xticks;
            xticklabels = out.xticklabels;
            On_task = out.order;
            
            Ar_rest = avg_restMat(On_task, On_task);
            %Ar_task = avg_corrMat(On_task, On_task);
        
        else    
            % use BCT function to order according to network affiliation
            [On_task Ar_task] = reorder_mod(avg_corrMat, parcellabels);
            % Reorder rest matrix in the same order as task
            Ar_rest = avg_restMat(On_task, On_task);

            %%
            % Map ordered parcels into scalars (for colormap)
            % Reorder parcels (to reflect the order on the correlogram
            networks_reordered = parcellabels(On_task);
            % Get the unique network association of each parcel
            uniqueparcels = unique(networks_reordered);
            % Map parcels to Xticks and Yticks
            mapObj = containers.Map(networks_reordered, 1:length(networks_reordered));
            % Get xtick and ytick labels
            [xticks index] = sort(cell2mat(mapObj.values));
            xticklabels = mapObj.keys;
            xticklabels = xticklabels(index); %sorted
        end
        
    end

    if strmatch(dtype, 'Volume')
        parcelfile = '/projects/ModalityControl/docs/scripts/Power11NewCoorOrder.txt';
        parcelimport = importdata(parcefile);
        parcellabels = parcelimport(1:264,5);
        

        % Make variables consistent in script
        Ar_rest = avg_restMat;
        %Ar_task = avg_corrMat;
    end 
end

% Outputs the correlation maps of the auditory versus visual task contrast
if step4 == 1
    modalConn = modalityCorrMaps(out_nogsr,subjNumStr, On_task, covariance);
    glmStats = GLMGroupAnalysis(out_nogsr, modalConn, On_task, taskRegression);
    if taskRegression==1 % most likely will not have to run this again, since I have run it once already and have retrived top 10 average activation ROIs in FPN (and DAN)
        overallTaskActivation_out = Cole2013SuppAnalysis(out_nogsr, cole2013_step);
    end
    cole2013_out = Cole2013SuppAnalysis(modalConn, cole2013_step)
end

    
% Plots data
if plot == 1
    figure
    Ar_task(1:n+1:n*n) = 0;
    imagesc(Ar_task)
    title('Task-Whole Timeseries')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7, 'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(Ar_task)),max(max((Ar_task)))))

    figure
    Ar_rest(1:n+1:n*n) = 0;
    imagesc(Ar_rest)
    title('Rest')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(Ar_rest)),max(max((Ar_rest)))))

    figure
    modalConn.corrMap_aud_avg(1:n+1:n*n) = 0;
    imagesc(modalConn.corrMap_aud_avg)
    title('Task - Auditory component only')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.corrMap_aud_avg)),max(max((modalConn.corrMap_aud_avg)))))

    figure
    modalConn.corrMap_vis_avg(1:n+1:n*n) = 0;
    imagesc(modalConn.corrMap_vis_avg)
    title('Task - Visual component only')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.corrMap_vis_avg)),max(max((modalConn.corrMap_vis_avg)))))
    
    figure
    modalConn.pval_Map(1:n+1:n*n) = 0;
    imagesc(modalConn.pval_Map)
    title('Auditory vs Visual p-value map')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.pval_Map)),max(max((modalConn.pval_Map)))))
   
    figure 
    modalConn.ttestBin_Map(1:n+1:n*n) = 0;
    imagesc(modalConn.ttestBin_Map)
    title('Binary T-Test Map')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.ttestBin_Map)),max(max((modalConn.ttestBin_Map)))))

    figure 
    modalConn.corrMap_diff(1:n+1:n*n) = 0;
    imagesc(modalConn.corrMap_diff)
    title('Visual - Auditory difference map')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.corrMap_diff)),max(max((modalConn.corrMap_diff)))))
    
    figure 
    modalConn.corrMap_diff_bin(1:n+1:n*n) = 0;
    imagesc(modalConn.corrMap_diff_bin)
    title('Visual - Auditory difference map; Thresholded')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colorbar()
    colormap(b2r(min(min(modalConn.corrMap_diff_bin)),max(max((modalConn.corrMap_diff_bin)))))

    figure
    modalConn.tstatPos(1:n+1:n*n) = 0;
    imagesc(modalConn.tstatPos)
    title('Positive Visual - Auditory T-stat')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colormap(b2r(min(min(modalConn.tstatPos)), max(max(modalConn.tstatPos))))
    colorbar()

    figure
    modalConn.tstatNeg(1:n+1:n*n) = 0;
    imagesc(modalConn.tstatNeg)
    title('Negative Visual - Auditory T-Stat')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize',7,'YDir', 'normal')
    set(gca, 'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colormap(b2r(min(min(modalConn.tstatNeg)), max(max(modalConn.tstatNeg))))
    colorbar()

    % Get correlations of significant ROIs
    indices = find(glmStats.ttest_bin);
    figure
    imagesc(modalConn.corrMap_diff(:,indices))
    tmp = parcellabels(On_task);
    xtick_tmp = tmp(indices);
    title('Whole-brain connectivity of significant task-evoked ROIs in Visual V Auditory responses')
    xlabel('Significant ROIs (labeled by network) with significant Visual v. Auditory task-evoked responses')
    ylabel('333 ROI regions sorted by network association')
    set(gca,'XTick', 1:length(indices),  'XTickLabel', xtick_tmp, 'YTick', xticks, 'YTickLabel', xticklabels, 'Fontsize', 7, 'YDir', 'normal')
    set(gca, 'xgrid', 'off', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k', 'LineWidth', 2)
    colormap(b2r(min(min(modalConn.corrMap_diff)), max(max(modalConn.corrMap_diff))))
    colorbar()

    
end


if writeSubjectDataToCSV == 1
    if strmatch(dtype, 'Surface')
        % First write GLM output to CSV
        outbasedir = '/projects/ModalityControl/data/';
        for i=1:length(subjNumStr)
            if excludeTaskEvokedResponseRegs==1
                optionStr = '_notaskModeled';
            else
                optionStr = '';
            end
            subj = subjNumStr{i};
            outdir = [outbasedir subj '/analysis/'];
            mkdir([outdir 'glmAnalysis'])
            restBetas = [outdir 'glmAnalysis/parcellated_restBetas.csv'];
            restResids = [outdir 'glmAnalysis/parcellated_restResids.csv'];
            taskBetas = [outdir 'glmAnalysis/parcellated_taskBetas' optionStr '.csv'];
            taskResids = [outdir 'glmAnalysis/parcellated_taskResids' optionStr '.csv'];
            binary_aud = [outdir 'glmAnalysis/binaryRegressors_aud.csv'];
            binary_vis = [outdir 'glmAnalysis/binaryRegressors_vis.csv'];
            rest_dtseries = [outdir 'glmAnalysis/rest_dtseries.csv'];
            task_dtseries = [outdir 'glmAnalysis/task_dtseries.csv'];

            csvwrite(restBetas, out_nogsr{i}.regionData_restBetas)
            csvwrite(restResids, out_nogsr{i}.regionData_restRegressed)
            csvwrite(taskBetas, out_nogsr{i}.regionData_taskBetas_byMod)
            csvwrite(taskResids, out_nogsr{i}.regionData_taskRegressed_byMod)
            csvwrite(binary_aud, out_nogsr{i}.binaryRegressors_aud)
            csvwrite(binary_vis, out_nogsr{i}.binaryRegressors_vis)
            csvwrite(rest_dtseries, out_nogsr{i}.dtseries_rest)
            csvwrite(task_dtseries, out_nogsr{i}.dtseries_task)
        end
    elseif strmatch(dtype, 'Volume')
        outbasedir = '/projects/ModalityControl/data/';
        for i=1:length(subjNumStr)
            subj = subjNumStr{i};
            outdir = [outbasedir subj '/analysis/'];
            mkdir([outdir 'glmAnalysis'])
            restBetas = [outdir 'glmAnalysis/parcellated_restBetas_PowerVol.csv'];
            restResids = [outdir 'glmAnalysis/parcellated_restResids_PowerVol.csv'];
            taskBetas = [outdir 'glmAnalysis/parcellated_taskBetas_PowerVol.csv'];
            taskResids = [outdir 'glmAnalysis/parcellated_taskResids_PowerVol.csv'];
            binary_aud = [outdir 'glmAnalysis/binaryRegressors_aud_PowerVol.csv'];
            binary_vis = [outdir 'glmAnalysis/binaryRegressors_vis_PowerVol.csv'];
            rest_dtseries = [outdir 'glmAnalysis/rest_dtseries_PowerVol.csv'];
            task_dtseries = [outdir 'glmAnalysis/task_dtseries_PowerVol.csv'];

            csvwrite(restBetas, out_nogsr{i}.regionData_restBetas)
            csvwrite(restResids, out_nogsr{i}.regionData_restRegressed)
            csvwrite(taskBetas, out_nogsr{i}.regionData_taskBetas_byMod)
            csvwrite(taskResids, out_nogsr{i}.regionData_taskRegressed_byMod)
            csvwrite(binary_aud, out_nogsr{i}.binaryRegressors_aud)
            csvwrite(binary_vis, out_nogsr{i}.binaryRegressors_vis)
            csvwrite(rest_dtseries, out_nogsr{i}.dtseries_rest)
            csvwrite(task_dtseries, out_nogsr{i}.dtseries_task)
        end
    end
end

if allVerticesGLM_subjectwise==1
    parfor (i=1:length(subjNumStr),7)
        disp(['Running on subject: ' subjNumStr{i}])
        tmp = surfaceGLM_forIPYNB(subjNumStr{i}, gsr);
    end
end

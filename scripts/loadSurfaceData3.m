function output = loadSurfaceData3(subj)
% Taku Ito
% 3/12/15
%
% This script loads in the surface data using the Gordon et al. 2014 333 ROI surface parcellation scheme
% 
% Parameters: subj ( must be in put with single quotations, i.e., as a string)

    % Get data directory based on the assumption that script is in projects/IndivRITL/docs/scripts/modalControl/
    basedir = ['/projects/IndivRITL/data/'];
    analysisdir = [basedir '/' subj '/analysis/'];
    datadir = [basedir subj '/MNINonLinear/Results'];
    mkdir([basedir subj]);
    mkdir(analysisdir);
    numTasks = 8;

    % Each dtseries is 333x581, so we will want to create an empty matrix of 333x4648 (8*581 = 4648)
    % We will create a 3d Matrix first, and reshape it into a 2d roi x dtseries after
    output_tmp = zeros(333,581,8);
    
    for task=1:numTasks

        % First need to parcellate each surface's Task${num}.dtseries.nii using workbench command
        inFile = [datadir '/Task' num2str(task) '/Task' num2str(task) '_Atlas.dtseries.nii'];
        labelFile = [basedir '/indivritl/Parcels/Parcels_LR.dlabel.nii'];
        outFile = [analysisdir '/Task' num2str(task) '_333Parcels.dtseries.nii'];
        % Run command through MATLAB's command line interface
        disp(['Running wb_command -cifti-parcellate to parcellate surface dtseries into 333 ROIs'])
        eval(['!wb_command -cifti-parcellate ' inFile ' ' labelFile ' COLUMN ' outFile])

        % Now, import the data to MATLAB (for a particular task) using ft_read_cifti in AnalysisTools
        disp(['Importing parcellated cifti file for task ' num2str(task)])
        data{task} = ft_read_cifti(outFile);

        % Now concatenate the time series across all 8 tasks. 
        % The variable data is a cell with 8 structs. To get the dtseries of each struct, we acceess them via data{task}.dtseries
        output_tmp(:,:,task) = data{task}.dtseries;
        % Demean each run
        disp(['Demeaning each run separately on task prior to concatenation...'])
        for roi=1:333
            output_tmp(roi,:,task) = output_tmp(roi,:,task) - mean(output_tmp(roi,:,task));
        end
    
    end

    % Parcellate and import rest data too
    inFile = [datadir '/Rest1/Rest1_Atlas.dtseries.nii'];
    labelFile = [basedir '/indivritl/Parcels/Parcels_LR.dlabel.nii'];
    outFile = [analysisdir '/Rest1_333Parcels.dtseries.nii'];

    % Run command through MATLAB's command line interface
    disp(['Running wb_command -cifti-parcellate to parcellate surface dtseries into 333 ROIs on Rest dataset'])
    eval(['!wb_command -cifti-parcellate ' inFile ' ' labelFile ' COLUMN ' outFile])

    % Now import the data to MATLAB using ft_read_cifti in AnalysisTools
    disp(['Importing parcellated cifti file for rest data'])
    rest = ft_read_cifti(outFile);

    % Reshape into two dimensions, 333 x 4648
    output.dtseries_task = reshape(output_tmp, [333,4648]);
    output.task = data;
    output.rest = rest;
    output.dtseries_rest = rest.dtseries;

end

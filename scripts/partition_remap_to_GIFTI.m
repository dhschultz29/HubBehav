function [template] = partition_remap_to_GIFTI(partitions,L_parcelPath,R_parcelPath,...
    num_parcels_per_hem,pathToTemplateGIFTI)

%Script author:
%Kaustubh Kulkarni
%
%Date: August 31, 2016

addpath('/projects/archive/GordonIndiv/docs/scripts/surface_parcellation/cifti-matlab-master/');
addpath('/projects/archive/GordonIndiv/docs/scripts/gifti-1.6/');

% Ensure that first half of partitions are left labels only
% and second half are right labels only
left_partitioned_labels = partitions(1:num_parcels_per_hem);
right_partitioned_labels = partitions(num_parcels_per_hem+1:2*num_parcels_per_hem);

% Load cifti labels
left_labels = ft_read_cifti(L_parcelPath);
right_labels = ft_read_cifti(R_parcelPath);

% Initialize gifti labels (32k vertices)
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

% Replace parcel label num with cluster num from louvain
for i=1:num_parcels_per_hem
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

% Rewrite to template gifti and save output
template = gifti(pathToTemplateGIFTI);
template.cdata = relabeled_left';
template.cdata = relabeled_right';
    
end
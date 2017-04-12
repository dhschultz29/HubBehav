%%%Calculate GVC for IndivRITL 64 task data%%%
%%%by Doug Schultz%%%
NetworkAssign=importdata(['/projects/AnalysisTools/netpartitions/ColeLabNetPartition_v1/parcel_network_assignments.txt']);
GVC_all_sub=zeros(360,100);
for subjNum=1:100
    connMat=All_subj_all_task_FC_r(:,:,:,subjNum);
    gvcVal=gvc(connMat);
    GVC_all_sub(:,subjNum)=gvcVal;
end

GVC_mean=nanmean(GVC_all_sub,2);

%%%GVC by network for each subj%%%
GVC_by_network=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=GVC_all_sub(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        GVC_by_network(netNum,subjNum)=NetGVC;
    end
end

%%map to gifti%%
%%import cifti files%%
left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = GVC_mean(1:180);
right_partitioned_labels = GVC_mean(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mean_GVC_R.func.gii','Base64Binary')
save(template,'mean_GVC_L.func.gii','Base64Binary')

%%%Calc GBC%%%

GBC_all_sub=zeros(360,100);
for subjNum=1:100
    FCmat=All_subj_all_rest_FC_partialr(:,:,subjNum);
    FCmat(FCmat==1)= nan;
    for region=1:360
        GBC_all_sub(region,subjNum)=nanmean(FCmat(region,:));
    end
end

GBC_all_sub_exclude_within_net=zeros(360,100);
for subjNum=1:100
    FCmat=All_subj_all_rest_FC_partialr(:,:,subjNum);
    FCmat(FCmat==1)= nan;
    for region=1:360
        communiteeNum=NetworkAssign(region,1);
        indexexcludewinetwork=(NetworkAssign~=communiteeNum);
        regionlvl=FCmat(region,:);
        excludewithinnetworkregions=nanmean(regionlvl(indexexcludewinetwork))';
        GBC_all_sub_exclude_within_net(region,subjNum)=excludewithinnetworkregions;
    end
end


    GBC_mean=nanmean(GBC_all_sub,2);
    GBC_mean_exclude_wi_net=nanmean(GBC_all_sub_exclude_within_net,2);
    
    GBC_by_network=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGBC=GBC_all_sub_exclude_within_net(:,subjNum);
        NetGBC=nanmean(subjGBC(NetworkAssign==netNum));
        GBC_by_network(netNum,subjNum)=NetGBC;
    end
end
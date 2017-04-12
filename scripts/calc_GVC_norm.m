%%%GVC Norm%%%

GVC_all_sub_norm=zeros(360,100);
for subjNum=1:100
    sublvl=GVC_all_sub(:,subjNum);
    for region=1:360
        regionlvl=sublvl(region,1);
        final=(regionlvl-nanmean(sublvl))/(nanstd(sublvl));
        GVC_all_sub_norm(region,subjNum)=final;
    end
end
   
GVC_by_network_norm=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=GVC_all_sub_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        GVC_by_network_norm(netNum,subjNum)=NetGVC;
    end
end

GVC_network_mean_norm=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_norm(netNum,:);
    netmean=nanmean(reglvl);
    GVC_network_mean_norm(1,netNum)=netmean;
end

GVC_network_std_norm=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_norm(netNum,:);
    netstd=nanstd(reglvl);
    GVC_network_std_norm(1,netNum)=netstd;
end

GVC_norm_corr_w_CPRO_R=zeros(360,1);
GVC_norm_corr_w_CPRO_P=zeros(360,1);
for regionNum=1:360
    data=GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_CPRO_bin),CPRO_Z(GVC_CPRO_bin));
    GVC_norm_corr_w_CPRO_R(regionNum,1)=R(1,2);
    GVC_norm_corr_w_CPRO_P(regionNum,1)=P(1,2);
end

GVC_norm_corr_w_Catell_R=zeros(360,1);
GVC_norm_corr_w_Catell_P=zeros(360,1);
for regionNum=1:360
    data=GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Catell_bin),Catell_Z(GVC_Catell_bin));
    GVC_norm_corr_w_Catell_R(regionNum,1)=R(1,2);
    GVC_norm_corr_w_Catell_P(regionNum,1)=P(1,2);
end

GVC_norm_corr_w_Raven_R=zeros(360,1);
GVC_norm_corr_w_Raven_P=zeros(360,1);
for regionNum=1:360
    data=GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Raven_bin),Raven_Z(GVC_Raven_bin));
    GVC_norm_corr_w_Raven_R(regionNum,1)=R(1,2);
    GVC_norm_corr_w_Raven_P(regionNum,1)=P(1,2);
end
    
left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = GVC_mean_norm(1:180);
right_partitioned_labels = GVC_mean_norm(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mean_GVC_norm_R.func.gii','Base64Binary')
save(template,'mean_GVC_norm_L.func.gii','Base64Binary')

%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = GVC_norm_corr_w_CPRO_R(1:180);
right_partitioned_labels = GVC_norm_corr_w_CPRO_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'GVC_norm_corr_w_CPRO_R.func.gii','Base64Binary')
save(template,'GVC_norm_corr_w_CPRO_L.func.gii','Base64Binary')
       

%%%
left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = GVC_norm_corr_w_Catell_R(1:180);
right_partitioned_labels = GVC_norm_corr_w_Catell_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'GVC_norm_corr_w_Catell_R.func.gii','Base64Binary')
save(template,'GVC_norm_corr_w_Catell_L.func.gii','Base64Binary')
       
%%%%%%
left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = GVC_norm_corr_w_Raven_R(1:180);
right_partitioned_labels = GVC_norm_corr_w_Raven_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'GVC_norm_corr_w_Raven_R.func.gii','Base64Binary')
save(template,'GVC_norm_corr_w_Raven_L.func.gii','Base64Binary')
 

GVC_by_network_norm_oldpart=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=GVC_all_sub_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign_old==netNum));
        GVC_by_network_norm_oldpart(netNum,subjNum)=NetGVC;
    end
end

GVC_network_mean_norm_oldpart=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_norm_oldpart(netNum,:);
    netmean=nanmean(reglvl);
    GVC_network_mean_norm_oldpart(1,netNum)=netmean;
end


GVC_network_std_norm_oldpart=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_norm_oldpart(netNum,:);
    netstd=nanstd(reglvl);
    GVC_network_std_norm_oldpart(1,netNum)=netstd;
end

%%%old GVC%%%newpartition%%%

GVC_by_network_newpart=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=GVC_all_sub(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        GVC_by_network_newpart(netNum,subjNum)=NetGVC;
    end
end

GVC_network_mean_newpart=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_newpart(netNum,:);
    netmean=nanmean(reglvl);
    GVC_network_mean_newpart(1,netNum)=netmean;
end


GVC_network_std_newpart=zeros(1,14);
for netNum=1:14
    reglvl=GVC_by_network_newpart(netNum,:);
    netstd=nanstd(reglvl);
    GVC_network_std_newpart(1,netNum)=netstd;
end

GVC_network_SEM_newpart=GVC_network_std_newpart/(sqrt(96));
All_subj_all_task_FC_norm=zeros(360,360,64,100);
normed=zeros(360,360,64);
for subjNum=1:100
    FCmat=All_subj_all_task_FC_r(:,:,:,subjNum);
    reshapedmat=reshape(FCmat,[8294400,1]);
    matmean=nanmean(reshapedmat);
    matstd=nanstd(reshapedmat);
    for taskNum=1:64
        for region1=1:360
            for region2=1:360
                normed(region1,region2,taskNum)=(FCmat(region1,region2,taskNum)-matmean)./matstd;
            end
        end
    end
    All_subj_all_task_FC_norm(:,:,:,subjNum)=normed;
end

    norm_GVC_all_sub=zeros(360,100);
for subjNum=1:100
    connMat=All_subj_all_task_FC_norm(:,:,:,subjNum);
    gvcVal=gvc(connMat);
    norm_GVC_all_sub(:,subjNum)=gvcVal;
end

norm_GVC_by_network=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=norm_GVC_all_sub(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        norm_GVC_by_network(netNum,subjNum)=NetGVC;
    end
end

norm_GVC_network_mean=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network(netNum,:);
    netmean=nanmean(reglvl);
    norm_GVC_network_mean(1,netNum)=netmean;
end


norm_GVC_network_std=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network(netNum,:);
    netstd=nanstd(reglvl);
    norm_GVC_network_std(1,netNum)=netstd;
end

norm_GVC_network_SEM=norm_GVC_network_std/(sqrt(99));

norm_GVC_corr_w_CPRO_R=zeros(360,1);
norm_GVC_corr_w_CPRO_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_CPRO_bin),CPRO_Z(GVC_CPRO_bin));
    norm_GVC_corr_w_CPRO_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_CPRO_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Catell_R=zeros(360,1);
norm_GVC_corr_w_Catell_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Catell_bin),Catell_Z(GVC_Catell_bin));
    norm_GVC_corr_w_Catell_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Catell_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Raven_R=zeros(360,1);
norm_GVC_corr_w_Raven_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Raven_bin),Raven_Z(GVC_Raven_bin));
    norm_GVC_corr_w_Raven_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Raven_P(regionNum,1)=P(1,2);
end

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_corr_w_CPRO_R(1:180);
right_partitioned_labels = norm_GVC_corr_w_CPRO_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_corr_CPRO_R.func.gii','Base64Binary')
save(template,'norm_GVC_corr_CPRO_L.func.gii','Base64Binary')

%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_corr_w_Catell_R(1:180);
right_partitioned_labels = norm_GVC_corr_w_Catell_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_corr_Catell_R.func.gii','Base64Binary')
save(template,'norm_GVC_corr_Catell_L.func.gii','Base64Binary')

%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_corr_w_Raven_R(1:180);
right_partitioned_labels = norm_GVC_corr_w_Raven_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_corr_Raven_R.func.gii','Base64Binary')
save(template,'norm_GVC_corr_Raven_L.func.gii','Base64Binary')

norm_GVC_all_sub_norm=zeros(360,100);
for subjNum=1:100
    sublvl=norm_GVC_all_sub(:,subjNum);
    for region=1:360
        regionlvl=sublvl(region,1);
        final=(regionlvl-nanmean(sublvl))/(nanstd(sublvl));
        norm_GVC_all_sub_norm(region,subjNum)=final;
    end
end
   
norm_GVC_by_network_norm=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=norm_GVC_all_sub_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        norm_GVC_by_network_norm(netNum,subjNum)=NetGVC;
    end
end

norm_GVC_network_mean_norm=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network_norm(netNum,:);
    netmean=nanmean(reglvl);
    norm_GVC_network_mean_norm(1,netNum)=netmean;
end

norm_GVC_network_std_norm=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network_norm(netNum,:);
    netstd=nanstd(reglvl);
    norm_GVC_network_std_norm(1,netNum)=netstd;
end

norm_GVC_network_SEM_norm=norm_GVC_network_std_norm/(sqrt(99));

norm_GVC_norm_corr_w_CPRO_R=zeros(360,1);
norm_GVC_norm_corr_w_CPRO_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_CPRO_bin),CPRO_Z(GVC_CPRO_bin));
    norm_GVC_norm_corr_w_CPRO_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_CPRO_P(regionNum,1)=P(1,2);
end

norm_GVC_norm_corr_w_Catell_R=zeros(360,1);
norm_GVC_norm_corr_w_Catell_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Catell_bin),Catell_Z(GVC_Catell_bin));
    norm_GVC_norm_corr_w_Catell_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_Catell_P(regionNum,1)=P(1,2);
end

norm_GVC_norm_corr_w_Raven_R=zeros(360,1);
norm_GVC_norm_corr_w_Raven_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Raven_bin),Raven_Z(GVC_Raven_bin));
    norm_GVC_norm_corr_w_Raven_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_Raven_P(regionNum,1)=P(1,2);
end

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_norm_corr_w_CPRO_R(1:180);
right_partitioned_labels = norm_GVC_norm_corr_w_CPRO_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_norm_corr_CPRO_R.func.gii','Base64Binary')
save(template,'norm_GVC_norm_corr_CPRO_L.func.gii','Base64Binary')

%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_norm_corr_w_Catell_R(1:180);
right_partitioned_labels = norm_GVC_norm_corr_w_Catell_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_norm_corr_Catell_R.func.gii','Base64Binary')
save(template,'norm_GVC_norm_corr_Catell_L.func.gii','Base64Binary')

%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_norm_corr_w_Raven_R(1:180);
right_partitioned_labels = norm_GVC_norm_corr_w_Raven_R(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'norm_GVC_norm_corr_Raven_R.func.gii','Base64Binary')
save(template,'norm_GVC_norm_corr_Raven_L.func.gii','Base64Binary')


%%%%Permutation test%%%%

Catell_GVC_R=zeros(360,10000);
Catell_GVC_P=zeros(360,10000);
for iteration=1:10000
iteration
rand=randperm(96);
Catell_shuff=Catell_Z_filt(rand);
behav=Catell_shuff;
for regNum=1:360
     tmp=norm_GVC_sub_norm_Catell(regNum,:);
     [R,P]=corrcoef(tmp,behav);
     Catell_GVC_R(regNum,iteration)=R(1,2);
     Catell_GVC_P(regNum,iteration)=P(1,2);
end
end

CPRO_GVC_R=zeros(360,10000);
CPRO_GVC_P=zeros(360,10000);
for iteration=1:10000
iteration
rand=randperm(94);
CPRO_shuff=CPRO_Z_filt(rand);
behav=CPRO_shuff;
for regNum=1:360
     tmp=norm_GVC_sub_norm_CPRO(regNum,:);
     [R,P]=corrcoef(tmp,behav);
     CPRO_GVC_R(regNum,iteration)=R(1,2);
     CPRO_GVC_P(regNum,iteration)=P(1,2);
end
end


Raven_GVC_R=zeros(360,10000);
Raven_GVC_P=zeros(360,10000);
for iteration=1:10000
iteration
rand=randperm(95);
Raven_shuff=Raven_Z_filt(rand);
behav=Raven_shuff;
for regNum=1:360
     tmp=norm_GVC_sub_norm_Raven(regNum,:);
     [R,P]=corrcoef(tmp,behav);
     Raven_GVC_R(regNum,iteration)=R(1,2);
     Raven_GVC_P(regNum,iteration)=P(1,2);
end

norm_GVC_norm_coor_w_Catell_adj_p=zeros(360,1);
norm_GVC_norm_coor_w_Catell_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=Catell_GVC_P(regNum,:);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Catell_thresh_p(regNum,1)=sorttempcellvec(500);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Catell_P(regNum,1)))./10000;
    norm_GVC_norm_coor_w_Catell_adj_p(regNum,1)=adjusted_p;
end

norm_GVC_norm_coor_w_CPRO_adj_p=zeros(360,1);
norm_GVC_norm_coor_w_CPRO_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=CPRO_GVC_P(regNum,:);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_CPRO_thresh_p(regNum,1)=sorttempcellvec(500);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_CPRO_P(regNum,1)))./10000;
    norm_GVC_norm_coor_w_CPRO_adj_p(regNum,1)=adjusted_p;
end

norm_GVC_norm_coor_w_Raven_adj_p=zeros(360,1);
norm_GVC_norm_coor_w_Raven_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=Raven_GVC_P(regNum,:);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Raven_thresh_p(regNum,1)=sorttempcellvec(500);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Raven_P(regNum,1)))./10000;
    norm_GVC_norm_coor_w_Raven_adj_p(regNum,1)=adjusted_p;
end

%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Catell(1:180);
right_partitioned_labels = mask_Catell(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Catell_R.func.gii','Base64Binary')
save(template,'mask_Catell_L.func.gii','Base64Binary')

%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_CPRO(1:180);
right_partitioned_labels = mask_CPRO(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_CPRO_R.func.gii','Base64Binary')
save(template,'mask_CPRO_L.func.gii','Base64Binary')

%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Raven(1:180);
right_partitioned_labels = mask_Raven(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Raven_R.func.gii','Base64Binary')
save(template,'mask_Raven_L.func.gii','Base64Binary')

%%%%%New permutation%%%%%%

%%reshuffle all FC%%
reshuff_GVC_Raven=zeros(360,10000,95);
for subjNum=1:95;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_GVC_singsubj_Raven=norm_GVC_sub_norm_Raven(:,subjNum);
        reshuff_GVC_Raven(:,iteration,subjNum)=norm_GVC_singsubj_Raven(randreg);
    end
end 

GVC_Raven_perm_R=zeros(360,10000);
GVC_Raven_perm_P=zeros(360,10000);
behav=Raven_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_GVC_Raven(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Raven_perm_R(regionNum,iteration)=R(1,2);
        GVC_Raven_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_Raven_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_Raven_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Raven_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Raven_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Raven_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_Raven_adj_p2(regNum,1)=adjusted_p;
end

mask_Raven2=norm_GVC_norm_corr_w_Raven_R.*(norm_GVC_norm_coor_w_Raven_adj_p2<0.05);

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Raven2(1:180);
right_partitioned_labels = mask_Raven2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Raven2_R.func.gii','Base64Binary')
save(template,'mask_Raven2_L.func.gii','Base64Binary')

%%%%%

reshuff_GVC_CPRO=zeros(360,10000,94);
for subjNum=1:94;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_GVC_singsubj_CPRO=norm_GVC_sub_norm_CPRO(:,subjNum);
        reshuff_GVC_CPRO(:,iteration,subjNum)=norm_GVC_singsubj_CPRO(randreg);
    end
end 

GVC_CPRO_perm_R=zeros(360,10000);
GVC_CPRO_perm_P=zeros(360,10000);
behav=CPRO_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_GVC_CPRO(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_CPRO_perm_R(regionNum,iteration)=R(1,2);
        GVC_CPRO_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_CPRO_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_CPRO_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_CPRO_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_CPRO_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_CPRO_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_CPRO_adj_p2(regNum,1)=adjusted_p;
end
mask_CPRO2=norm_GVC_norm_corr_w_CPRO_R.*(norm_GVC_norm_coor_w_CPRO_adj_p2<0.05);

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_CPRO2(1:180);
right_partitioned_labels = mask_CPRO2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_CPRO2_R.func.gii','Base64Binary')
save(template,'mask_CPRO2_L.func.gii','Base64Binary')

%%%%%

reshuff_GVC_Catell=zeros(360,10000,96);
for subjNum=1:96;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_GVC_singsubj_Catell=norm_GVC_sub_norm_Catell(:,subjNum);
        reshuff_GVC_Catell(:,iteration,subjNum)=norm_GVC_singsubj_Catell(randreg);
    end
end 

GVC_Catell_perm_R=zeros(360,10000);
GVC_Catell_perm_P=zeros(360,10000);
behav=Catell_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_GVC_Catell(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Catell_perm_R(regionNum,iteration)=R(1,2);
        GVC_Catell_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_Catell_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_Catell_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Catell_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Catell_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Catell_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_Catell_adj_p2(regNum,1)=adjusted_p;
end
mask_Catell2=norm_GVC_norm_corr_w_Catell_R.*(norm_GVC_norm_coor_w_Catell_adj_p2<0.05);

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Catell2(1:180);
right_partitioned_labels = mask_Catell2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Catell2_R.func.gii','Base64Binary')
save(template,'mask_Catell2_L.func.gii','Base64Binary')


%%%%Modularity%%%
Ci_all=zeros(360,100);
Q_all=zeros(100,1);
for subjNum=1:100
    connMat=All_subj_all_rest_FC_GSR_r(:,:,subjNum);
    connMat(isnan(connMat))=1;
    connMat_thresh=threshold_proportional(connMat,0.5);
    [Ci Q]=modularity_und(connMat_thresh,1);
    Ci_all(:,subjNum)=Ci(:,1);
    Q_all(subjNum,1)=Q(1,1);
end

upperTriangle2=(tril(ones(size(upperTriangle)))==0);
mean_rest_FC=zeros(100,1);
for subjNum=1:100
    connMat=All_subj_all_rest_FC_GSR_r(:,:,subjNum);
    connMat1D=connMat(upperTriangle2);
    FCmean=nanmean(connMat1D);
    mean_rest_FC(subjNum,1)=FCmean;
end

%%%%corr w gen factor%%%%
norm_GVC_norm_corr_w_Gen_factor_R=zeros(360,1);
norm_GVC_norm_corr_w_Gen_factor_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Gen_factor_bin),Gen_factor_filt);
    norm_GVC_norm_corr_w_Gen_factor_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_Gen_factor_P(regionNum,1)=P(1,2);
end

reshuff_GVC_Gen_factor=zeros(360,10000,99);
for subjNum=1:99;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_GVC_singsubj_Gen_factor=norm_GVC_sub_norm_Gen_factor(:,subjNum);
        reshuff_GVC_Gen_factor(:,iteration,subjNum)=norm_GVC_singsubj_Gen_factor(randreg);
    end
end 

GVC_Gen_factor_perm_R=zeros(360,10000);
GVC_Gen_factor_perm_P=zeros(360,10000);
behav=Gen_factor_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_GVC_Gen_factor(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Gen_factor_perm_R(regionNum,iteration)=R(1,2);
        GVC_Gen_factor_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_Gen_factor_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_Gen_factor_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Gen_factor_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Gen_factor_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Gen_factor_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_Gen_factor_adj_p2(regNum,1)=adjusted_p;
end

mask_Gen_factor2=norm_GVC_norm_corr_w_Gen_factor_R.*(norm_GVC_norm_coor_w_Gen_factor_adj_p2<0.05);


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Gen_factor2(1:180);
right_partitioned_labels = mask_Gen_factor2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Gen_factor2_R.func.gii','Base64Binary')
save(template,'mask_Gen_factor2_L.func.gii','Base64Binary')

%%%%%%%%%%%%%%

norm_GVC_norm_corr_w_DCS_R=zeros(360,1);
norm_GVC_norm_corr_w_DCS_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_DCS_bin),DCS_Z_filt);
    norm_GVC_norm_corr_w_DCS_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_DCS_P(regionNum,1)=P(1,2);
end

reshuff_GVC_DCS=zeros(360,10000,91);
for subjNum=1:91;
    for iteration=1:10000;
        (iteration*subjNum)/(91*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_DCS=norm_GVC_sub_norm_DCS(:,subjNum);
        reshuff_GVC_DCS(:,iteration,subjNum)=norm_GVC_singsubj_DCS(randreg);
    end
end 

GVC_DCS_perm_R=zeros(360,10000);
GVC_DCS_perm_P=zeros(360,10000);
behav=DCS_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (iteration*regionNum)/3600000
        GVCsamp=squeeze(reshuff_GVC_DCS(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_DCS_perm_R(regionNum,iteration)=R(1,2);
        GVC_DCS_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_DCS_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_DCS_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_DCS_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_DCS_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_DCS_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_DCS_adj_p2(regNum,1)=adjusted_p;
end

mask_DCS2=norm_GVC_norm_corr_w_DCS_R.*(norm_GVC_norm_coor_w_DCS_adj_p2<0.05);


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_DCS2(1:180);
right_partitioned_labels = mask_DCS2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_DCS2_R.func.gii','Base64Binary')
save(template,'mask_DCS2_L.func.gii','Base64Binary')

%%%%%%%%%%%%%%%%%

norm_GVC_norm_corr_w_Flanker_R=zeros(360,1);
norm_GVC_norm_corr_w_Flanker_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Flanker_bin),Flanker_Z_filt);
    norm_GVC_norm_corr_w_Flanker_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_Flanker_P(regionNum,1)=P(1,2);
end

reshuff_GVC_Flanker=zeros(360,10000,92);
for subjNum=1:92;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(92*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_Flanker=norm_GVC_sub_norm_Flanker(:,subjNum);
        reshuff_GVC_Flanker(:,iteration,subjNum)=norm_GVC_singsubj_Flanker(randreg);
    end
end 

GVC_Flanker_perm_R=zeros(360,10000);
GVC_Flanker_perm_P=zeros(360,10000);
behav=Flanker_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_Flanker(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Flanker_perm_R(regionNum,iteration)=R(1,2);
        GVC_Flanker_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_norm_coor_w_Flanker_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_Flanker_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Flanker_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Flanker_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Flanker_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_Flanker_adj_p2(regNum,1)=adjusted_p;
end

mask_Flanker2=norm_GVC_norm_corr_w_Flanker_R.*(norm_GVC_norm_coor_w_Flanker_adj_p2<0.05);


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Flanker2(1:180);
right_partitioned_labels = mask_Flanker2(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Flanker2_R.func.gii','Base64Binary')
save(template,'mask_Flanker2_L.func.gii','Base64Binary')

%%%%Calculate a matrix of FC STD values%%%

All_subj_all_task_FC_STD=zeros(360,360,100);
for subjNum=1:100
    FCmat=All_subj_all_task_FC_norm(:,:,:,subjNum);
    for reg1=1:360
        for reg2=1:360
            All_subj_all_task_FC_STD(reg1,reg2,subjNum)=std(FCmat(reg1,reg2,:));
        end
    end
end

%%%Now calculate correlation with gen_factor at each connection%%%
STD_Corr_w_Gen_factor_R=zeros(360,360);
STD_Corr_w_Gen_factor_P=zeros(360,360);
STDmat=All_subj_all_task_FC_STD(:,:,(GVC_bin));
for reg1=1:360
    for reg2=1:360
        (((reg1-1)*360)+reg2)/(360*360)
        [R,P]=corrcoef(Gen_factor_filt,STDmat(reg1,reg2,:));
        STD_Corr_w_Gen_factor_R(reg1,reg2)=R(1,2);
        STD_Corr_w_Gen_factor_P(reg1,reg2)=P(1,2);
    end
end


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = FPN2_region_FC_STD_corr_w_Gen_factor(1:180);
right_partitioned_labels = FPN2_region_FC_STD_corr_w_Gen_factor(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'FPN2_region_FC_STD_corr_w_Gen_factor_R.func.gii','Base64Binary')
save(template,'FPN2_region_FC_STD_corr_w_Gen_factor_L.func.gii','Base64Binary')


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = SM_region_FC_STD_corr_w_Gen_factor(1:180);
right_partitioned_labels = SM_region_FC_STD_corr_w_Gen_factor(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'SM2_region_FC_STD_corr_w_Gen_factor_R.func.gii','Base64Binary')
save(template,'SM2_region_FC_STD_corr_w_Gen_factor_L.func.gii','Base64Binary')

norm_GVC_norm_by_network=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=norm_GVC_all_sub_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        norm_GVC_norm_by_network(netNum,subjNum)=NetGVC;
    end
end

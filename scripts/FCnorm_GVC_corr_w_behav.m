norm_GVC_corr_w_CPRO_R=zeros(360,1);
norm_GVC_corr_w_CPRO_P=zeros(360,1);
for regionNum=1:360
    (regionNum)/2160
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_CPRO_bin),CPRO_Z(GVC_CPRO_bin));
    norm_GVC_corr_w_CPRO_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_CPRO_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Catell_R=zeros(360,1);
norm_GVC_corr_w_Catell_P=zeros(360,1);
for regionNum=1:360
    (regionNum+360)/2160
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Catell_bin),Catell_Z(GVC_Catell_bin));
    norm_GVC_corr_w_Catell_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Catell_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Raven_R=zeros(360,1);
norm_GVC_corr_w_Raven_P=zeros(360,1);
for regionNum=1:360
    (regionNum+720)/2160
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Raven_bin),Raven_Z(GVC_Raven_bin));
    norm_GVC_corr_w_Raven_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Raven_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Gen_factor_R=zeros(360,1);
norm_GVC_corr_w_Gen_factor_P=zeros(360,1);
for regionNum=1:360
    (regionNum+1080)/2160
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Gen_factor_bin),Gen_factor_filt);
    norm_GVC_corr_w_Gen_factor_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Gen_factor_P(regionNum,1)=P(1,2);
end
norm_GVC_corr_w_DCS_R=zeros(360,1);
norm_GVC_corr_w_DCS_P=zeros(360,1);
for regionNum=1:360
    (regionNum+1440)/2160
    data=norm_GVC_all_sub(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_DCS_bin),DCS_Z_filt);
    norm_GVC_corr_w_DCS_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_DCS_P(regionNum,1)=P(1,2);
end

norm_GVC_corr_w_Flanker_R=zeros(360,1);
norm_GVC_corr_w_Flanker_P=zeros(360,1);
for regionNum=1:360
    (regionNum+1800)/2160
    data=norm_GVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Flanker_bin),Flanker_Z_filt);
    norm_GVC_corr_w_Flanker_R(regionNum,1)=R(1,2);
    norm_GVC_corr_w_Flanker_P(regionNum,1)=P(1,2);
end

%%%%%%%%%%%%%%%%clean%%%%%%%%%%%%%%%%%%%%%%%

norm_GVC_sub_Flanker=norm_GVC_all_sub(:,GVC_Flanker_bin==1);

norm_GVC_sub_Gen_factor=norm_GVC_all_sub(:,GVC_Gen_factor_bin==1);

norm_GVC_sub_Catell=norm_GVC_all_sub(:,GVC_Catell_bin==1);

norm_GVC_sub_CPRO=norm_GVC_all_sub(:,GVC_CPRO_bin==1);

norm_GVC_sub_DCS=norm_GVC_all_sub(:,GVC_DCS_bin==1);

norm_GVC_sub_Raven=norm_GVC_all_sub(:,GVC_Raven_bin==1);

%%%%%%%%%%%%%%%%Permutations%%%%%%%%%%%%%%%%

reshuff_GVC_Flanker=zeros(360,10000,92);
for subjNum=1:92;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(92*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_Flanker=norm_GVC_sub_Flanker(:,subjNum);
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

norm_GVC_coor_w_Flanker_adj_p=zeros(360,1);
norm_GVC_coor_w_Flanker_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Flanker_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_Flanker_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_Flanker_P(regNum,1)))./3600000;
    norm_GVC_coor_w_Flanker_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_Flanker=norm_GVC_corr_w_Flanker_R.*(norm_GVC_coor_w_Flanker_adj_p<0.05);

%%%%

reshuff_GVC_Catell=zeros(360,10000,96);
for subjNum=1:96;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(96*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_Catell=norm_GVC_sub_Catell(:,subjNum);
        reshuff_GVC_Catell(:,iteration,subjNum)=norm_GVC_singsubj_Catell(randreg);
    end
end 

GVC_Catell_perm_R=zeros(360,10000);
GVC_Catell_perm_P=zeros(360,10000);
behav=Catell_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_Catell(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Catell_perm_R(regionNum,iteration)=R(1,2);
        GVC_Catell_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_coor_w_Catell_adj_p=zeros(360,1);
norm_GVC_coor_w_Catell_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Catell_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_Catell_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_Catell_P(regNum,1)))./3600000;
    norm_GVC_coor_w_Catell_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_Catell=norm_GVC_corr_w_Catell_R.*(norm_GVC_coor_w_Catell_adj_p<0.05);

%%%

reshuff_GVC_CPRO=zeros(360,10000,94);
for subjNum=1:94;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(94*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_CPRO=norm_GVC_sub_CPRO(:,subjNum);
        reshuff_GVC_CPRO(:,iteration,subjNum)=norm_GVC_singsubj_CPRO(randreg);
    end
end 

GVC_CPRO_perm_R=zeros(360,10000);
GVC_CPRO_perm_P=zeros(360,10000);
behav=CPRO_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_CPRO(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_CPRO_perm_R(regionNum,iteration)=R(1,2);
        GVC_CPRO_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_coor_w_CPRO_adj_p=zeros(360,1);
norm_GVC_coor_w_CPRO_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_CPRO_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_CPRO_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_CPRO_P(regNum,1)))./3600000;
    norm_GVC_coor_w_CPRO_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_CPRO=norm_GVC_corr_w_CPRO_R.*(norm_GVC_coor_w_CPRO_adj_p<0.05);

%%%

reshuff_GVC_DCS=zeros(360,10000,91);
for subjNum=1:91;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(91*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_DCS=norm_GVC_sub_Flanker(:,subjNum);
        reshuff_GVC_DCS(:,iteration,subjNum)=norm_GVC_singsubj_DCS(randreg);
    end
end 

GVC_DCS_perm_R=zeros(360,10000);
GVC_DCS_perm_P=zeros(360,10000);
behav=DCS_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_DCS(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_DCS_perm_R(regionNum,iteration)=R(1,2);
        GVC_DCS_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_coor_w_DCS_adj_p=zeros(360,1);
norm_GVC_coor_w_DCS_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_DCS_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_DCS_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_DCS_P(regNum,1)))./3600000;
    norm_GVC_coor_w_DCS_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_DCS=norm_GVC_corr_w_DCS_R.*(norm_GVC_coor_w_DCS_adj_p<0.05);

%%%

reshuff_GVC_Gen_factor=zeros(360,10000,99);
for subjNum=1:99;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(99*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_Gen_factor=norm_GVC_sub_Gen_factor(:,subjNum);
        reshuff_GVC_Gen_factor(:,iteration,subjNum)=norm_GVC_singsubj_Gen_factor(randreg);
    end
end 

GVC_Gen_factor_perm_R=zeros(360,10000);
GVC_Gen_factor_perm_P=zeros(360,10000);
behav=Gen_factor_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_Gen_factor(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Gen_factor_perm_R(regionNum,iteration)=R(1,2);
        GVC_Gen_factor_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_coor_w_Gen_factor_adj_p=zeros(360,1);
norm_GVC_coor_w_Gen_factor_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Gen_factor_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_Gen_factor_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_Gen_factor_P(regNum,1)))./3600000;
    norm_GVC_coor_w_Gen_factor_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_Gen_factor=norm_GVC_corr_w_Gen_factor_R.*(norm_GVC_coor_w_Gen_factor_adj_p<0.05);

%%%

reshuff_GVC_Raven=zeros(360,10000,95);
for subjNum=1:95;
    for iteration=1:10000;
        ((subjNum*10000)+iteration)/(95*10000)
        randreg=randperm(360);
        norm_GVC_singsubj_Raven=norm_GVC_sub_Raven(:,subjNum);
        reshuff_GVC_Raven(:,iteration,subjNum)=norm_GVC_singsubj_Raven(randreg);
    end
end 

GVC_Raven_perm_R=zeros(360,10000);
GVC_Raven_perm_P=zeros(360,10000);
behav=Raven_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        (((regionNum-1)*10000)+iteration)/3600000
        GVCsamp=squeeze(reshuff_GVC_Raven(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Raven_perm_R(regionNum,iteration)=R(1,2);
        GVC_Raven_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_GVC_coor_w_Raven_adj_p=zeros(360,1);
norm_GVC_coor_w_Raven_thresh_p=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Raven_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_coor_w_Raven_thresh_p(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_corr_w_Raven_P(regNum,1)))./3600000;
    norm_GVC_coor_w_Raven_adj_p(regNum,1)=adjusted_p;
end

mask_FCnorm_Raven=norm_GVC_corr_w_Raven_R.*(norm_GVC_coor_w_Raven_adj_p<0.05);

%%%%%%%Plot%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_Flanker(181:360);
right_partitioned_labels = mask_FCnorm_Flanker(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_Flanker_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_Flanker_L.func.gii','Base64Binary')

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_Catell(181:360);
right_partitioned_labels = mask_FCnorm_Catell(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_Catell_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_Catell_L.func.gii','Base64Binary')

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_CPRO(181:360);
right_partitioned_labels = mask_FCnorm_CPRO(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_CPRO_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_CPRO_L.func.gii','Base64Binary')

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_DCS(181:360);
right_partitioned_labels = mask_FCnorm_DCS(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_DCS_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_DCS_L.func.gii','Base64Binary')

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_Gen_factor(181:360);
right_partitioned_labels = mask_FCnorm_Gen_factor(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_Gen_factor_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_Gen_factor_L.func.gii','Base64Binary')

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_FCnorm_Raven(181:360);
right_partitioned_labels = mask_FCnorm_Raven(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_FCnorm_Raven_R.func.gii','Base64Binary')
save(template,'mask_FCnorm_Raven_L.func.gii','Base64Binary')


left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = norm_GVC_norm_mean(181:360);
right_partitioned_labels = norm_GVC_norm_mean(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_norm_GVC_norm_mean_R.func.gii','Base64Binary')
save(template,'mask_norm_GVC_norm_mean_L.func.gii','Base64Binary')

%%%%Task to rest comparison%%%

upperTriangle=(tril(ones(360,360))==0);

Rest_to_Task_sim=zeros(100,64);
for subjNum=1:100
for taskNum=1:64
	Resttmp=All_subj_all_rest_FC_Z(:,:,subjNum);
	Tasktmp=All_subj_all_task_FC_Z(:,:,taskNum,subjNum);
	[R,P]=corrcoef(Resttmp(upperTriangle),Tasktmp(upperTriangle))
	Rest_to_Task_sim(subjNum,taskNum)=R(1,2);
end
end

Task_to_avetask_sim=zeros(100,64);
for subjNum=1:100
    for taskNum=1:64
        Tasktmp=All_subj_all_task_FC_Z(:,:,taskNum,subjNum);
        AveTasktmp=nanmean(All_subj_all_task_FC_Z(:,:,NewtaskNum~=taskNum,subjNum),3);
	    [R,P]=corrcoef(Tasktmp(upperTriangle),AveTasktmp(upperTriangle))
	    Task_to_avetask_sim(subjNum,taskNum)=R(1,2);
    end
end

All_subj_all_mean_task_FC_Z=nanmean(All_subj_all_task_FC_Z,3);

Rest_to_aveTask_sim=zeros(100,1);
for subjNum=1:100
	Resttmp=All_subj_all_rest_FC_Z(:,:,subjNum);
	Tasktmp=All_subj_all_mean_task_FC_Z(:,:,subjNum);
	[R,P]=corrcoef(Resttmp(upperTriangle),Tasktmp(upperTriangle))
	Rest_to_aveTask_sim(subjNum,1)=R(1,2);
end

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = TaskbetaT_FDRp_05(181:360);
right_partitioned_labels = TaskbetaT_FDRp_05(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'TaskbetaT_FDRp_05_R.func.gii','Base64Binary')
save(template,'TaskbetaT_FDRp_05_L.func.gii','Base64Binary')

%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = TaskbetaT_FDRp_005(181:360);
right_partitioned_labels = TaskbetaT_FDRp_005(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'TaskbetaT_FDRp_005_R.func.gii','Base64Binary')
save(template,'TaskbetaT_FDRp_005_L.func.gii','Base64Binary')
    norm_outGVC_all_sub=zeros(360,100);
    
    for regNum=1:360
        NetworkAssigntmp=NetworkAssign;
        NetworkAssigntmp(regNum,1)=29;
        %tmpConnmat=zeros(size(NetworkAssigntmp~=NetworkAssign(regNum,1),NetworkAssigntmp~=NetworkAssign(regNum,1),64,100));
        tmpConnmat=All_subj_all_task_FC_norm(NetworkAssigntmp~=NetworkAssign(regNum,1),NetworkAssigntmp~=NetworkAssign(regNum,1),:,:);
        regind=find(NetworkAssigntmp~=NetworkAssign(regNum,1));
        for subjNum=1:100
            PercentDone=((((regNum-1)*100)+subjNum)/36000)*100
            connMat=tmpConnmat(:,:,:,subjNum);
            gvcVal=gvc(connMat);
            for regind1=1:size(regind);
            norm_outGVC_all_sub(regind(regind1),subjNum)=gvcVal((regind1),1);
            end
        end
        clear NetworkAssigntmp
        clear tmpConnmat
        clear connMat
        clear regind
    end
            
    
norm_outGVC_all_sub_norm=zeros(360,100);
for subjNum=1:100
    sublvl=norm_outGVC_all_sub(:,subjNum);
    for region=1:360
        regionlvl=sublvl(region,1);
        final=(regionlvl-nanmean(sublvl))/(nanstd(sublvl));
        norm_outGVC_all_sub_norm(region,subjNum)=final;
    end
end

norm_outGVC_by_network_norm=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=norm_outGVC_all_sub_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        norm_outGVC_by_network_norm(netNum,subjNum)=NetGVC;
    end
end

norm_outGVC_network_mean_norm=zeros(1,14);
for netNum=1:14
    reglvl=norm_outGVC_by_network_norm(netNum,:);
    netmean=nanmean(reglvl);
    norm_outGVC_network_mean_norm(1,netNum)=netmean;
end

norm_outGVC_network_std=zeros(1,14);
for netNum=1:14
    reglvl=norm_outGVC_by_network_norm(netNum,:);
    netstd=nanstd(reglvl);
    norm_outGVC_network_std(1,netNum)=netstd;
end

norm_outGVC_network_SEM=norm_outGVC_network_std/(sqrt(99));

norm_outGVC_norm_corr_w_Gen_factor_R=zeros(360,1);
norm_outGVC_norm_corr_w_Gen_factor_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_Gen_factor_bin),Gen_factor_filt);
    norm_outGVC_norm_corr_w_Gen_factor_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_Gen_factor_P(regionNum,1)=P(1,2);
end

norm_outGVC_norm_corr_w_Catell_R=zeros(360,1);
norm_outGVC_norm_corr_w_Catell_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_Catell_bin),Catell_Z_filt);
    norm_outGVC_norm_corr_w_Catell_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_Catell_P(regionNum,1)=P(1,2);
end

norm_outGVC_norm_corr_w_CPRO_R=zeros(360,1);
norm_outGVC_norm_corr_w_CPRO_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_CPRO_bin),CPRO_Z_filt);
    norm_outGVC_norm_corr_w_CPRO_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_CPRO_P(regionNum,1)=P(1,2);
end

norm_outGVC_norm_corr_w_DCS_R=zeros(360,1);
norm_outGVC_norm_corr_w_DCS_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_DCS_bin),DCS_Z_filt);
    norm_outGVC_norm_corr_w_DCS_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_DCS_P(regionNum,1)=P(1,2);
end

norm_outGVC_norm_corr_w_Flanker_R=zeros(360,1);
norm_outGVC_norm_corr_w_Flanker_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_Flanker_bin),Flanker_Z_filt);
    norm_outGVC_norm_corr_w_Flanker_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_Flanker_P(regionNum,1)=P(1,2);
end

norm_outGVC_norm_corr_w_Raven_R=zeros(360,1);
norm_outGVC_norm_corr_w_Raven_P=zeros(360,1);
for regionNum=1:360
    data=norm_outGVC_all_sub_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(outGVC_Raven_bin),Raven_Z_filt);
    norm_outGVC_norm_corr_w_Raven_R(regionNum,1)=R(1,2);
    norm_outGVC_norm_corr_w_Raven_P(regionNum,1)=P(1,2);
end

reshuff_outGVC_Gen_factor=zeros(360,10000,99);
for subjNum=1:99;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_Gen_factor=norm_outGVC_sub_norm_Gen_factor(:,subjNum);
        reshuff_outGVC_Gen_factor(:,iteration,subjNum)=norm_outGVC_singsubj_Gen_factor(randreg);
    end
end 

reshuff_outGVC_Catell=zeros(360,10000,96);
for subjNum=1:96;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_Catell=norm_outGVC_sub_norm_Catell(:,subjNum);
        reshuff_outGVC_Catell(:,iteration,subjNum)=norm_outGVC_singsubj_Catell(randreg);
    end
end 

reshuff_outGVC_CPRO=zeros(360,10000,94);
for subjNum=1:94;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_CPRO=norm_outGVC_sub_norm_CPRO(:,subjNum);
        reshuff_outGVC_CPRO(:,iteration,subjNum)=norm_outGVC_singsubj_CPRO(randreg);
    end
end 

reshuff_outGVC_DCS=zeros(360,10000,91);
for subjNum=1:91;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_DCS=norm_outGVC_sub_norm_DCS(:,subjNum);
        reshuff_outGVC_DCS(:,iteration,subjNum)=norm_outGVC_singsubj_DCS(randreg);
    end
end 

reshuff_outGVC_Flanker=zeros(360,10000,92);
for subjNum=1:92;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_Flanker=norm_outGVC_sub_norm_Flanker(:,subjNum);
        reshuff_outGVC_Flanker(:,iteration,subjNum)=norm_outGVC_singsubj_Flanker(randreg);
    end
end 

reshuff_outGVC_Raven=zeros(360,10000,95);
for subjNum=1:95;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_outGVC_singsubj_Raven=norm_outGVC_sub_norm_Raven(:,subjNum);
        reshuff_outGVC_Raven(:,iteration,subjNum)=norm_outGVC_singsubj_Raven(randreg);
    end
end 

outGVC_Gen_factor_perm_R=zeros(360,10000);
outGVC_Gen_factor_perm_P=zeros(360,10000);
behav=Gen_factor_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_Gen_factor(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_Gen_factor_perm_R(regionNum,iteration)=R(1,2);
        outGVC_Gen_factor_perm_P(regionNum,iteration)=P(1,2);
    end
end

outGVC_Catell_perm_R=zeros(360,10000);
outGVC_Catell_perm_P=zeros(360,10000);
behav=Catell_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_Catell(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_Catell_perm_R(regionNum,iteration)=R(1,2);
        outGVC_Catell_perm_P(regionNum,iteration)=P(1,2);
    end
end

outGVC_CPRO_perm_R=zeros(360,10000);
outGVC_CPRO_perm_P=zeros(360,10000);
behav=CPRO_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_CPRO(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_CPRO_perm_R(regionNum,iteration)=R(1,2);
        outGVC_CPRO_perm_P(regionNum,iteration)=P(1,2);
    end
end

outGVC_DCS_perm_R=zeros(360,10000);
outGVC_DCS_perm_P=zeros(360,10000);
behav=DCS_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_DCS(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_DCS_perm_R(regionNum,iteration)=R(1,2);
        outGVC_DCS_perm_P(regionNum,iteration)=P(1,2);
    end
end

outGVC_Flanker_perm_R=zeros(360,10000);
outGVC_Flanker_perm_P=zeros(360,10000);
behav=Flanker_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_Flanker(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_Flanker_perm_R(regionNum,iteration)=R(1,2);
        outGVC_Flanker_perm_P(regionNum,iteration)=P(1,2);
    end
end

outGVC_Raven_perm_R=zeros(360,10000);
outGVC_Raven_perm_P=zeros(360,10000);
behav=Raven_Z_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_outGVC_Raven(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        outGVC_Raven_perm_R(regionNum,iteration)=R(1,2);
        outGVC_Raven_perm_P(regionNum,iteration)=P(1,2);
    end
end

norm_outGVC_norm_coor_w_Gen_factor_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_Gen_factor_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_Gen_factor_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_Gen_factor_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_Gen_factor_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_Gen_factor_adj_p2(regNum,1)=adjusted_p;
end
mask_Gen_factor_outGVC=norm_outGVC_norm_corr_w_Gen_factor_R.*(norm_outGVC_norm_coor_w_Gen_factor_adj_p2<0.05);

norm_outGVC_norm_coor_w_Catell_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_Catell_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_Catell_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_Catell_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_Catell_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_Catell_adj_p2(regNum,1)=adjusted_p;
end
mask_Catell_outGVC=norm_outGVC_norm_corr_w_Catell_R.*(norm_outGVC_norm_coor_w_Catell_adj_p2<0.05);

norm_outGVC_norm_coor_w_CPRO_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_CPRO_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_CPRO_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_CPRO_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_CPRO_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_CPRO_adj_p2(regNum,1)=adjusted_p;
end
mask_CPRO_outGVC=norm_outGVC_norm_corr_w_CPRO_R.*(norm_outGVC_norm_coor_w_CPRO_adj_p2<0.05);

norm_outGVC_norm_coor_w_DCS_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_DCS_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_DCS_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_DCS_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_DCS_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_DCS_adj_p2(regNum,1)=adjusted_p;
end
mask_DCS_outGVC=norm_outGVC_norm_corr_w_DCS_R.*(norm_outGVC_norm_coor_w_DCS_adj_p2<0.05);

norm_outGVC_norm_coor_w_Flanker_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_Flanker_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_Flanker_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_Flanker_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_Flanker_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_Flanker_adj_p2(regNum,1)=adjusted_p;
end
mask_Flanker_outGVC=norm_outGVC_norm_corr_w_Flanker_R.*(norm_outGVC_norm_coor_w_Flanker_adj_p2<0.05);

norm_outGVC_norm_coor_w_Raven_adj_p2=zeros(360,1);
norm_outGVC_norm_coor_w_Raven_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(outGVC_Raven_perm_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_outGVC_norm_coor_w_Raven_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_outGVC_norm_corr_w_Raven_P(regNum,1)))./3600000;
    norm_outGVC_norm_coor_w_Raven_adj_p2(regNum,1)=adjusted_p;
end
mask_Raven_outGVC=norm_outGVC_norm_corr_w_Raven_R.*(norm_outGVC_norm_coor_w_Raven_adj_p2<0.05);

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Gen_factor_outGVC(181:360);
right_partitioned_labels = mask_Gen_factor_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Gen_factor_outGVC_R.func.gii','Base64Binary')
save(template,'mask_Gen_factor_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Catell_outGVC(181:360);
right_partitioned_labels = mask_Catell_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Catell_outGVC_R.func.gii','Base64Binary')
save(template,'mask_Catell_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_CPRO_outGVC(181:360);
right_partitioned_labels = mask_CPRO_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_CPRO_outGVC_R.func.gii','Base64Binary')
save(template,'mask_CPRO_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_DCS_outGVC(181:360);
right_partitioned_labels = mask_DCS_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_DCS_outGVC_R.func.gii','Base64Binary')
save(template,'mask_DCS_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Flanker_outGVC(181:360);
right_partitioned_labels = mask_Flanker_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Flanker_outGVC_R.func.gii','Base64Binary')
save(template,'mask_Flanker_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Raven_outGVC(181:360);
right_partitioned_labels = mask_Raven_outGVC(1:180);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Raven_outGVC_R.func.gii','Base64Binary')
save(template,'mask_Raven_outGVC_L.func.gii','Base64Binary')

%%%%%%%%%%%%%
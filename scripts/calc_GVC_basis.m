All_subj_all_task_FC_basis_norm=zeros(360,360,64,100);
normed=zeros(360,360,64);
for subjNum=1:100
    FCmat=All_subj_all_task_FC_basis_r(:,:,:,subjNum);
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
    All_subj_all_task_FC_basis_norm(:,:,:,subjNum)=normed;
end

    norm_GVC_all_sub_basis=zeros(360,100);
for subjNum=1:100
    connMat=All_subj_all_task_FC_basis_norm(:,:,:,subjNum);
    gvcVal=gvc(connMat);
    norm_GVC_all_sub_basis(:,subjNum)=gvcVal;
end

norm_GVC_all_sub_basis_norm=zeros(360,100);
for subjNum=1:100
    sublvl=norm_GVC_all_sub_basis(:,subjNum);
    for region=1:360
        regionlvl=sublvl(region,1);
        final=(regionlvl-nanmean(sublvl))/(nanstd(sublvl));
        norm_GVC_all_sub_basis_norm(region,subjNum)=final;
    end
end

norm_GVC_by_network_basis_norm=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGVC=norm_GVC_all_sub_basis_norm(:,subjNum);
        NetGVC=nanmean(subjGVC(NetworkAssign==netNum));
        norm_GVC_by_network_basis_norm(netNum,subjNum)=NetGVC;
    end
end

norm_GVC_network_mean_basis_norm=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network_basis_norm(netNum,:);
    netmean=nanmean(reglvl);
    norm_GVC_network_mean_basis_norm(1,netNum)=netmean;
end

norm_GVC_network_std_basis=zeros(1,14);
for netNum=1:14
    reglvl=norm_GVC_by_network_basis_norm(netNum,:);
    netstd=nanstd(reglvl);
    norm_GVC_network_std_basis(1,netNum)=netstd;
end

norm_GVC_network_SEM_basis=norm_GVC_network_std_basis/(sqrt(99));

norm_GVC_norm_corr_w_Gen_factor_basis_R=zeros(360,1);
norm_GVC_norm_corr_w_Gen_factor_basis_P=zeros(360,1);
for regionNum=1:360
    data=norm_GVC_all_sub_basis_norm(regionNum,:);
    data=data';
    [R,P]=corrcoef(data(GVC_Gen_factor_bin),Gen_factor_filt);
    norm_GVC_norm_corr_w_Gen_factor_basis_R(regionNum,1)=R(1,2);
    norm_GVC_norm_corr_w_Gen_factor_basis_P(regionNum,1)=P(1,2);
end

reshuff_GVC_Gen_factor_basis=zeros(360,10000,99);
for subjNum=1:99;
    for iteration=1:10000;
        iteration
        randreg=randperm(360);
        norm_GVC_singsubj_Gen_factor_basis=norm_GVC_sub_norm_Gen_factor_basis(:,subjNum);
        reshuff_GVC_Gen_factor_basis(:,iteration,subjNum)=norm_GVC_singsubj_Gen_factor_basis(randreg);
    end
end 

%%%%%%%%%%%%%%%%%%

GVC_Gen_factor_perm_basis_R=zeros(360,10000);
GVC_Gen_factor_perm_basis_P=zeros(360,10000);
behav=Gen_factor_filt;
for regionNum=1:360;
    for iteration=1:10000;
        GVCsamp=squeeze(reshuff_GVC_Gen_factor_basis(regionNum,iteration,:));
        [R,P]=corrcoef(GVCsamp,behav);
        GVC_Gen_factor_perm_basis_R(regionNum,iteration)=R(1,2);
        GVC_Gen_factor_perm_basis_P(regionNum,iteration)=P(1,2);
    end
end

%%%%%%%%%%%%%%%%%%

norm_GVC_norm_coor_w_Gen_factor_basis_adj_p2=zeros(360,1);
norm_GVC_norm_coor_w_Gen_factor_basis_thresh_p2=zeros(360,1);
for regNum=1:360
    tempcellvec=reshape(GVC_Gen_factor_perm_basis_P(:,:),[3600000,1]);
    sorttempcellvec=sort(tempcellvec);
    norm_GVC_norm_coor_w_Gen_factor_basis_thresh_p2(regNum,1)=sorttempcellvec(180000);
    adjusted_p=(sum(tempcellvec<norm_GVC_norm_corr_w_Gen_factor_basis_P(regNum,1)))./3600000;
    norm_GVC_norm_coor_w_Gen_factor_basis_adj_p2(regNum,1)=adjusted_p;
end

mask_Gen_factor2_basis=norm_GVC_norm_corr_w_Gen_factor_basis_R.*(norm_GVC_norm_coor_w_Gen_factor_basis_adj_p2<0.05);

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mask_Gen_factor2_basis(1:180);
right_partitioned_labels = mask_Gen_factor2_basis(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mask_Gen_factor2_basis_R.func.gii','Base64Binary')
save(template,'mask_Gen_factor2_basis_L.func.gii','Base64Binary')

%%%%%%%%%%%%%

left_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
right_labels = ft_read_cifti(['/mnt/engram/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']);
relabeled_left = zeros(1,length(left_labels.data));
relabeled_right = zeros(1,length(right_labels.data));

left_partitioned_labels = mean_norm_GVC_norm_basis(1:180);
right_partitioned_labels = mean_norm_GVC_norm_basis(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'mean_norm_GVC_norm_basis_R.func.gii','Base64Binary')
save(template,'mean_norm_GVC_norm_basis_L.func.gii','Base64Binary')



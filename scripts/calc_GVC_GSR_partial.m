%%%Calculate GVC for IndivRITL 64 task data%%%
%%%by Doug Schultz%%%
NetworkAssign=importdata(['/projects/AnalysisTools/netpartitions/ColeLabNetPartition_v1/parcel_network_assignments.txt']);
GVC_all_sub=zeros(360,100);
for subjNum=1:100
    connMat=atanh(All_subj_all_task_FC_GSR_r(:,:,:,subjNum));
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

left_partitioned_labels = CESD_GBC_corr_parcel_r(1:180);
right_partitioned_labels = CESD_GBC_corr_parcel_r(181:360);

for i=1:180
    relabeled_left(logical(left_labels.data==i)) = left_partitioned_labels(i);
    relabeled_right(logical(right_labels.data==i)) = right_partitioned_labels(i);    
end

template = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_L_0.6.func.gii']);
template2 = gifti(['/projects/archive/GordonIndiv/data/results/parcellation_output_140824/final_parcellation_R_0.6.func.gii']);
template.cdata = relabeled_left';
template2.cdata = relabeled_right';
save(template2,'CESD_GBC_corr_R.func.gii','Base64Binary')
save(template,'CESD_GBC_corr_L.func.gii','Base64Binary')

%%%Calc GBC%%%

GBC_all_sub=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_partialr(:,:,subjNum));
    FCmat(FCmat==Inf)= nan;
    for region=1:360
        GBC_all_sub(region,subjNum)=nanmean(FCmat(region,:));
    end
end

GBC_all_sub_exclude_within_net_=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_partialr(:,:,subjNum));
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

RestFC_GSR_corr_w_CESD_R=zeros(360,360);
RestFC_GSR_corr_w_CESD_P=zeros(360,360);
for i=1:360
    for j=1:360
         tmp=squeeze(atanh(All_subj_all_rest_FC_GSR_r(i,j,:)));
         behav=CESD_boxcox_Z;
         [R,P]=corrcoef(tmp(subselect),behav(subselect));
         RestFC_GSR_corr_w_CESD_R(i,j)=R(1,2);
         RestFC_GSR_corr_w_CESD_P(i,j)=P(1,2);
     end
end 

GBC_by_network_resids=zeros(14,100);
X=zeros(100,1);
for net=1:14
    X=GBC_by_network(net,:)';
    Y=net;
    allnet=1:14;
    regressors=GBC_by_network(allnet~=Y,:)';
    reg_resid=regstats(X, regressors, 'linear', {'r','beta','rsquare'});
    GBC_by_network_resids(net,:)=reg_resid.r'
end

CESD_Ttest_T=zeros(360,360);
CESD_Ttest_P=zeros(360,360);
for regionA=1:360
        for regionB=1:360
                [h,p,ci,stats]=ttest(atanh(All_subj_all_rest_FC_GSR_r(regionA,regionB,CESD_High)),atanh(All_subj_all_rest_FC_GSR_r(regionA,regionB,CESD_Low)));
                CESD_Ttest_T(regionA,regionB)=stats.tstat;
                CESD_Ttest_P(regionA,regionB)=p;
            end
end
        
CESD_GBC_corr_parcel_r=zeros(360,1);
CESD_GBC_corr_parcel_p=zeros(360,1);

for i=1:360
         tmp=squeeze(GBC_all_sub_exclude_within_net(i,CESD_mask));
         behav=CESD_all_behav(CESD_mask);
         [R,P]=corrcoef(tmp,behav);
         CESD_GBC_corr_parcel_r(i,:)=R(1,2);
         CESD_GBC_corr_parcel_p(i,:)=P(1,2);
end   

GBC_all_sub_within_net=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==Inf)= nan;
    for region=1:360
        communiteeNum=NetworkAssign(region,1);
        indexexcludewinetwork=(NetworkAssign==communiteeNum);
        regionlvl=FCmat(region,:);
        excludewithinnetworkregions=nanmean(regionlvl(indexexcludewinetwork))';
        GBC_all_sub_within_net(region,subjNum)=excludewithinnetworkregions;
    end
end

GBC_wi_network=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGBC=GBC_all_sub_within_net(:,subjNum);
        NetGBC=nanmean(subjGBC(NetworkAssign==netNum));
        GBC_wi_network(netNum,subjNum)=NetGBC;
    end
end

%%positive only GBC%%%



GBC_all_sub_exclude_within_net_pos=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==Inf)= nan;
    for ROI=1:360
        for ROI2=1:360
    FCmat_pos_mask(ROI,ROI2) = FCmat(ROI,ROI2)>0;
        end
    end
    for region=1:360
        communiteeNum=NetworkAssign(region,1);
        indexexcludewinetwork=(NetworkAssign==communiteeNum);
        indexexcludewinetwork2=indexexcludewinetwork(FCmat_pos_mask(region,:));
        regionlvl=FCmat(region,:);
        regionlvl2=regionlvl(FCmat_pos_mask(region,:));
        excludewithinnetworkregions=nanmean(regionlvl2(indexexcludewinetwork2))';
        GBC_all_sub_exclude_within_net_pos(region,subjNum)=excludewithinnetworkregions;
    end
end

GBC_mean_exclude_wi_network_pos=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGBC=GBC_all_sub_exclude_within_net_pos(:,subjNum);
        NetGBC=nanmean(subjGBC(NetworkAssign==netNum));
        GBC_mean_exclude_wi_network_pos(netNum,subjNum)=NetGBC;
    end
end

meanGBC_by_net_exclude_pos=zeros(14,1);
for net=1:14
meanGBC_by_net_exclude_pos(net,1)=nanmean(GBC_mean_exclude_wi_network_pos(net,:),2);
end

%%%GBC neg only%%%


GBC_all_sub_exclude_within_net_neg=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==Inf)= nan;
    for ROI=1:360
        for ROI2=1:360
    FCmat_neg_mask(ROI,ROI2) = FCmat(ROI,ROI2)<0;
        end
    end
    for region=1:360
        communiteeNum=NetworkAssign(region,1);
        indexexcludewinetwork=(NetworkAssign==communiteeNum);
        indexexcludewinetwork2=indexexcludewinetwork(FCmat_neg_mask(region,:));
        regionlvl=FCmat(region,:);
        regionlvl2=regionlvl(FCmat_neg_mask(region,:));
        excludewithinnetworkregions=nanmean(regionlvl2(indexexcludewinetwork2))';
        GBC_all_sub_exclude_within_net_neg(region,subjNum)=excludewithinnetworkregions;
    end
end

GBC_mean_exclude_wi_network_neg=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGBC=GBC_all_sub_exclude_within_net_neg(:,subjNum);
        NetGBC=nanmean(subjGBC(NetworkAssign==netNum));
        GBC_mean_exclude_wi_network_neg(netNum,subjNum)=NetGBC;
    end
end

meanGBC_by_net_exclude_neg=zeros(14,1);
for net=1:14
meanGBC_by_net_exclude_neg(net,1)=nanmean(GBC_mean_exclude_wi_network_neg(net,:),2);
end

GBC_all_sub_exclude_within_net_DMN_FPN_together=zeros(360,100);
for subjNum=1:100
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==1)= nan;
    for region=1:360
        communiteeNum=NetworkAssign(region,1);
        indexexcludewinetwork=(NetworkAssign~=communiteeNum);
        regionlvl=FCmat(region,:);
        excludewithinnetworkregions=nanmean(regionlvl(indexexcludewinetwork2))';
        GBC_all_sub_exclude_within_net_DMN_FPN_together(region,subjNum)=excludewithinnetworkregions;
    end
end


    GBC_mean=nanmean(GBC_all_sub,2);
    GBC_mean_exclude_wi_net_DMN_FPN_together=nanmean(GBC_all_sub_exclude_within_net_DMN_FPN_together,2);
    
    GBC_by_network_DMN_FPN_together=zeros(14,100);
for subjNum=1:100
    for netNum=1:14
        subjGBC=GBC_all_sub_exclude_within_net_DMN_FPN_together(:,subjNum);
        NetGBC=nanmean(subjGBC(NetworkAssign==netNum));
        GBC_by_network_DMN_FPN_together(netNum,subjNum)=NetGBC;
    end
end

FPN_DMN_FC_mean=zeros(100,1);

for subjNum=1:100;
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==1)= nan;
    FPN_DMN_only=FCmat(NetworkAssign==6,NetworkAssign==7);
    FPN_DMN_only2=reshape(FPN_DMN_only,[],1);
    FPN_DMN_FC_mean(subjNum,1)=nanmean(FPN_DMN_only2);
end
    
%%GBC_net_by_net%%

GBC_net_by_net_mean=zeros(100,196);



for subjNum=1:100;
    comb=1
    FCmat=atanh(All_subj_all_rest_FC_GSR_r(:,:,subjNum));
    FCmat(FCmat==Inf)= nan;
        for netnum1=1:14
            for netnum2=1:14
    sample=FCmat(NetworkAssign==netnum1,NetworkAssign==netnum2);
    sample2=reshape(sample,[],1);
    newmean(1,comb)=nanmean(sample2);
    GBC_net_by_net_mean(subjNum,comb)=nanmean(sample2);
    comb=comb+1
            end
        end
    clear sample;
    clear sample2;
end


%%%Calculate FC for 64 different tasks (& rest) from IndivRITL preprocessed
%%%data
%%%script by Doug Schultz
%%%takes "output" from Preproc_PostHCPMinimalPreproc_v1_2_DHS as input

All_subj_all_task_FC_GSR_r=zeros(360,360,64,100);
All_subj_all_task_FC_GSR_p=zeros(360,360,64,100);

for subjNum=1:100
    subjNum
    all_task_FC=zeros(360,360,64);
        %%make scrubbed task timing regressors
        scrubbedTaskTiming=output.taskTiming{subjNum}.taskdesignmat_hrf(logical(output.taskGLMVars{subjNum}.temporalmask),:);
        scrubbedTaskTiming_binary=scrubbedTaskTiming>0;
        scrubbedTaskTiming_binary_trans=scrubbedTaskTiming_binary';
    for taskNum=1:64
        taskKeep=scrubbedTaskTiming_binary_trans(taskNum,:);
        for regNum=1:360
        tasktimepoints(regNum,:)=output.taskGLMVars{subjNum}.fMRI_resids(regNum,:);
        mask=tasktimepoints(regNum,:);
        tasktimepoints_bytask(regNum,:)=mask(taskKeep);
        clear tasktimepoints
        clear mask
        end
        [TaskFC_corr_r,TaskFC_corr_p]=corrcoef(tasktimepoints_bytask');
        All_subj_all_task_FC_GSR_r(:,:,taskNum,subjNum)=TaskFC_corr_r;
        All_subj_all_task_FC_GSR_p(:,:,taskNum,subjNum)=TaskFC_corr_p;
        clear tasktimepoints_bytask
    end
end
       
All_subj_all_rest_FC_GSR_r=zeros(360,360,100);
All_subj_all_rest_FC_GSR_p=zeros(360,360,100);
        
for subjNum=1:100
    subjNum
        [RestFC_corr_r,RestFC_corr_p]=corrcoef(output.restGLMVars{subjNum}.fMRI_resids');
        All_subj_all_rest_FC_GSR_r(:,:,subjNum)=RestFC_corr_r;
        All_subj_all_rest_FC_GSR_p(:,:,subjNum)=RestFC_corr_p;
end        
        
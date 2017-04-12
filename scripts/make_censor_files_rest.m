%script to make censor files
cd /projects/IndivRITL/data/motionfiles
subjNums='013 014 016 017 018 021 023 024 025 026 027 028 030 031 032 033 034 035 037 038 039 040 042 043 045 046 047 048 049 050 052 053 055 056 057 058 062 063 066 067 068 069 070 072 074 075 076 077 078 079 081 085 086 087 088 090 092 093 094 095 097 098 099 101 102 104 105 106 108 109 110 111 112 114 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 134 137 138 139 140 141';
%subjNums='013';
subjNumStr = strread(subjNums, '%s', 'delimiter', ' ');
runstart=load(['runstart_rest.1D']);

for subjNum=1:length(subjNumStr);
    filename=subjNumStr{subjNum,1};
    data=load ([num2str(filename) '_Rest_FD.1D']);
    temp=data<0.5;
    temp2=double(temp);
    scrubbed_task=ones(1070,1);
    indexx=find(temp2(:,1)==0);
    indexx1=(indexx-1);
    indexx1(indexx1==0)=1;
    indexx2=(indexx+1);
    indexx3=(indexx+2);
    for ind=1:length(indexx)
        scrubbed_task(indexx(ind,1),1)=0;
    end
    for ind1=1:length(indexx1)
        scrubbed_task(indexx1(ind1,1),1)=0;
    end
    for ind2=1:length(indexx2)
        scrubbed_task(indexx2(ind2,1),1)=0;
    end
     for ind3=1:length(indexx3)
        scrubbed_task(indexx3(ind3,1),1)=0;
     end
    scrubbed_task_final=(scrubbed_task(1:1070,:));
    scrubbed_task_final2=(scrubbed_task_final.*(runstart));
    dlmwrite(([num2str(filename) '_Rest_censor.1D']),scrubbed_task_final2);
    filename
    percent_censored=1-(mean(scrubbed_task_final2(:,1)));
    percent_censored
end
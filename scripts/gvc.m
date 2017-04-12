function gvcVal = gvc(connMat)

%Function to calculate global variability coefficient, based on Cole M.W., Reynolds J.R., Power J.D., Repovs G., Anticevic A., Braver T.S. (2013). "Multi-task connectivity reveals flexible hubs for adaptive task control". Nature Neuroscience; 2013 Sep;16(9):1348-55. doi:10.1038/nn.3470.
%
%Input: NxNxS connectivity matrix, where N is the number of nodes and S is the number of states of those nodes.
%Output: The global variability coefficient for each node. The coefficient is calculated as the average of the standard deviations of a given node's connection strengths across all other nodes.
%
%By Michael W. Cole, mwcole@mwcole.net
%Version 0.1

n=size(connMat,1);                  %number of nodes
I=eye(n)~=0;                        %logical identity matrix
Imult=repmat(I,[1 1 size(connMat,3)]);    %full logical identity matrix
connMat(Imult)=NaN;

gvcVal=nanmean(std(connMat,0,3),2);
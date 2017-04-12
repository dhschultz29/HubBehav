numsubjects=12;
%input the number of subjects
nummeasures=20;
%input the number of measures
examplemat=randi(10,numsubjects,nummeasures);
%create a sample matrix of random integers (in this case 12 subjects by 20
%measures)
corrtable_r=zeros(nummeasures,nummeasures);
%create an empty table for us to fill with correlation values
corrtable_p=zeros(nummeasures,nummeasures);
%create an empty table for us to fill with p values
for firstmeasure=1:nummeasures
    for secondmeasure=1:nummeasures
        [R,P]=corrcoef(examplemat(:,firstmeasure),examplemat(:,secondmeasure));
        corrtable_r(firstmeasure,secondmeasure)=R(1,2);
        corrtable_r(secondmeasure,firstmeasure)=R(1,2);
        corrtable_p(firstmeasure,secondmeasure)=P(1,2);
        corrtable_p(secondmeasure,firstmeasure)=P(1,2);
    end
end
%run the correlations looping through the measures and output the results
%to our previously blank table



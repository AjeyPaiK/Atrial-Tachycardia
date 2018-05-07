                      %Source Code%
%---------------------------------------------------------------%
% Matlab code for the diagnosis of Atrial Tachycardia.          %
% Authors    :  Ajey Pai K, Chaitra Suresh, Arpith B.           %
% Affiliation:  Dept of E&CE, NMAMIT, Nitte, Karnataka, India.  %
% Outsourced by Piitech systems Pvt ltd, Bengaluru,India.       %
%---------------------------------------------------------------%
tic
clear;
clc;
load PwaveDB.mat;
load STPdb.mat;
load ATPdb.mat;
%-----------------------------------------------------------------------------------%
% Converting a record to wfdb format for pre-processing 
%-----------------------------------------------------------------------------------%
Fs=122;         %122Hz for Records provided. 360Hz for MITDB records.
adu='mV';
info='testrec'; %The output header file, data file will be generated in this name.
M=csvread('AF1mV-Lead2.csv');
mat2wfdb(M,info,Fs,[],adu,info);
%-----------------------------------------------------------------------------------%
% Pre-Processing using the WFDB commands to extract various parameters
%-----------------------------------------------------------------------------------%
sqrs(info,1220,1,[],[]);
ecgpuwave(info,'test',1,1220,'qrs',[],[]);
h=tach(info,'qrs',1220,1);
HR=mean(h);
[wavesstart,type,subtype,channel,nums,comments]=rdann(info,'test',[],[],[],'(');
[wavesend,type,subtype,channel,nume,comments]=rdann(info,'test',[],[],[],')');
[tm,signal]=rdsamp(info,[],1000);
%-----------------------------------------------------------------------------------%
%Initialising variables
%-----------------------------------------------------------------------------------%
p=0;                %to store samples at which p peaks occur
rs=0;              %to store samples at which r peaks occur
PPinterval=0;       %to store P-P interval
PRinterval=0;       %to store P-R interval
count1=0;
count2=0;
count3=0;
%-----------------------------------------------------------------------------------%
% Extracting sample values at which P peaks and R peaks occur
%------------------------------------------------------------%
for i=1:length(nums)
    if nums(i)==0
       count1=count1+1;
       p(count1)=wavesstart(i);
    end
    if nums(i)==1
        count2=count2+1;
        rs(count2)=wavesstart(i);
    end
end
%------------------------------------------------------------%
%Finding P-P intervals and storing them in vector PPinterval 
%------------------------------------------------------------%
for i=1:2:length(p)
    if i~=length(p)
        PPinterval(i)=(p(i+1)-p(i))/122;
    end
end
PPinterval=PPinterval(:,1:2:end);
%------------------------------------------------------------------%
%Finding the Atrial Rate (For Atrial Tachycardia Atrial Rate>100BPM)
%------------------------------------------------------------------%
ATR= mean(PPinterval);
ATR=60/ATR;
%------------------------------------------------%
%Equalising the lengths of vector p and vector r 
%------------------------------------------------%
if length(p)>length(rs)
    rs=[rs zeros(1,length(p)-length(rs))];
end
if length(rs)>length(p)
    p=[p zeros(1,length(rs)-length(p))];
end
%------------------------------------------------%
%generating vector pr to find P-R interval
%------------------------------------------------%
pr=zeros(1,2*length(p));
pr(1,1:2:end)=p;
pr(1,2:2:end)=rs;
%------------------------------------------------------------%
% Finding P-R intervals and storing them in vector PRinterval
%------------------------------------------------------------%
for i=1:2:length(pr)
    if i~=length(pr)
        PRinterval(i)=(pr(i+1)-pr(i))/122;
    end
end
PRinterval=PRinterval(:,1:2:end);
PRinterval=mean(PRinterval);
%-----------------------------------------------------------------------%
%Finding QRS Duration (For Atrial Tachycardia, QRS duration< 0.08s-0.12s
%-----------------------------------------------------------------------%
for i=1:length(nume)
    if nume(i)==1
        count3=count3+1;
        re(count3)=wavesend(i);
    end
end
rs=[rs;zeros(1,length(rs))];
rs=rs(:);
rs=rs';
rs(1,2:2:end)=re;
clear re;
for i=1:2:length(rs)
    if i~=length(rs)
        QRSinterval(i)=(rs(i+1)-rs(i))/122;
    end
end
QRSinterval=QRSinterval(:,1:2:end);
QRSinterval=mean(QRSinterval);
%------------------------------------------------------------%
% Extracting P-waves from the input record for Correlation   %
%------------------------------------------------------------%
count3=0;
count4=0;
for i=1:length(wavesstart)
    if nums(i)==0
       count3=count3+1;
       pwave(count3)=wavesstart(i);
    end
end
for i=1:length(nums)
    if nums(i)==0
       count4=count4+1;
       pe(count4)=wavesend(i);
    end
end
pwave=[pwave;zeros(1,length(pwave))];
pwave=pwave(:);
pwave=pwave';
pwave(1,2:2:end)=pe;
clear pe;
%--------------------------------------------------------------------%
% Cross Correlation of p-waves with different types of P-waves
%
% 1. Cross Correlation of p-waves with Standard Normal sinus P-waves
%--------------------------------------------------------------------%
featuresofp(50,50)=0;
d=0;
for i=1:2:length(pwave)
    if pwave(i)<=500
        d=signal(pwave(i):pwave(i+1));
        featuresofp(1:length(d),i)=d;
    else
        break;
    end
    clear d;
end
featuresofp(:,all(~any(featuresofp),1))=[];
[a,b]=size(featuresofp);
[e,f]=size(normalpdb);
normalpdb=normalpdb';
normalpdb=[normalpdb zeros(1,length(featuresofp)-length(normalpdb))];
normalpdb=normalpdb';
for i=1:b
    for j=1:f
        R=corrcoef(normalpdb(:,j),featuresofp(:,i));
        NPcorrv(i)=mean(min(R));
    end
    NPcorr=mean(NPcorrv);
end
%---------------------------------------------------------------%
% 2. Cross Correlation of p-waves with Sinus Tachycardia p-waves
%---------------------------------------------------------------%
[e,f]=size(STPdb);
STPdb=STPdb';
STPdb=[STPdb zeros(1,length(featuresofp)-length(STPdb))];
STPdb=STPdb';
for i=1:b
    for j=1:f
        R=corrcoef(STPdb(:,j),featuresofp(:,i));
        STPcorrv(i)=mean(min(R));
    end
    STPcorr=mean(STPcorrv);
end
%---------------------------------------------------------------%
% 3. Cross Correlation of p-waves with Atrial Tachycardia P-waves
%---------------------------------------------------------------%
[e,f]=size(ATPdb);
ATPdb=ATPdb';
ATPdb=[ATPdb zeros(1,length(featuresofp)-length(ATPdb))];
ATPdb=ATPdb';
for i=1:b
    for j=1:f
        R=corrcoef(ATPdb(:,j),featuresofp(:,i));
        ATPcorrv(i)=mean(min(R));
    end
    ATPcorr=mean(ATPcorrv);
end
clc;
%------------------------------------------------------------------%
%Printing Results
%------------------------------------------------------------------%
fprintf('The Heart Rate is:%.2f BPM\n\n',HR);
fprintf('The Atrial Rate is:%.2f BPM\n\n',ATR);
fprintf('The average QRS duration is:%.2f s\n\n',QRSinterval);
fprintf('Percentage Correlation with Normal Sinus P-waves is: %.2f Percent\n\n',NPcorr*100);
fprintf('Percentage Correlation with Atrial Tachycardia P-waves is: %.2f Percent\n\n',ATPcorr*100);
fprintf('Percentage Correlation with Sinus Tachycardia P-waves is: %.2f Percent\n\n',STPcorr*100);
%------------------------------------------------------------------%
%Printing Decision
%------------------------------------------------------------------%
if (HR<100&&ATR<100)
    fprintf('Atrial Tachycardia not present1\n\n');
    
elseif (ATR>120)
    if QRSinterval>0.12
        fprintf('Atrial Tachycardia with Delta wave is present2\n\n');
    else
        fprintf('Atrial Tachycardia is present3\n\n');
    end
    
elseif(HR>110&&ATR>100&&abs(STPcorr*100)>(abs(NPcorr*100)))
        fprintf('Atrial Tachycardia is not Present4\n\n');
elseif(HR>110&&ATR>100&&abs(ATPcorr*100)>(abs(STPcorr*100))&&abs(ATPcorr*100)>=abs(NPcorr*100))
        fprintf('Atrial Tachycardia is Present5\n\n');
elseif(HR>100&&ATR<100)&&(abs(ATPcorr*100)>abs(STPcorr*100))&&(abs(ATPcorr*100)>abs(NPcorr*100))
    if QRSinterval>0.12
        fprintf('Atrial Tachycardia with Delta wave is present6\n\n');
    else
        fprintf('Atrial Tachycardia is present7\n\n');
    end
else
    fprintf('Atrial Tachycardia is not present8\n\n');
end
toc
%End of Code%
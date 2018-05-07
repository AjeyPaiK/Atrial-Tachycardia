figure(3);
[a,b]=size(featuresofp);
t=0:0.01:0.22;
for i=1:b
    plot(t,featuresofp(1:23,i),'r');
    hold on;
end
title('P-wave features of an Atrial Tachycardia Record');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
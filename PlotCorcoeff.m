ATPcorrv=abs(ATPcorrv).*100;
STPcorrv=abs(STPcorrv).*100;
NPcorrv=abs(NPcorrv).*100;
plot(ATPcorrv,'g-*');
hold on
plot(STPcorrv,'b-o');
plot(NPcorrv,'r-x');
hold off
xlabel('sample number (n)');
ylabel('Correlation Coefficient (%)');
legend('with ATPdb','with STPdb','with NSRPdb','Location','east');

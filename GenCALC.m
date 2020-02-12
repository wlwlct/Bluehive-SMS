%%
function [CALC,TOTSIG]=GenCALC(tau,comparerange)

expvalue=1;
TOTSIG=0;
CALC=zeros(comparerange+1,1);
for i=1:comparerange+1
    %i basically could be understand as bin number,so there would be total
    %4000 ps range
t= 8*(i-1);%8 means bin size to be 8 ps.
CALC(i,1)=expvalue*exp(-t/tau);
TOTSIG=TOTSIG+CALC(i,1);%running sum to calc area

end

end
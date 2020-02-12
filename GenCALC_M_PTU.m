function [CALC,TOTSIG]=GenCALC_M(tau,comparerange_CALC)

expvalue=1;
TOTSIG=0;
t=0:comparerange_CALC;
    %i basically could be understand as bin number,so there would be total
    %4000 ps range
CALC=transpose(expvalue*exp(-t./tau));
TOTSIG=sum(CALC);

end
function DISTORED=GenDISTORED_M(CALC,IRFI_sub,TOTSIG,comparerange_CALC)
A=0;
DSUM=0;

DISTORED=zeros(comparerange_CALC,1);
for ii=1:comparerange_CALC
    DISTORED(ii,1)=sum(flipud(CALC(1:ii,1)).*IRFI_sub(1:ii,1));
DSUM=DSUM+DISTORED(ii,1);
    

end

DISTORED=DISTORED.*(TOTSIG/DSUM);
[DISTORED_max,~]=max(DISTORED);
DISTORED=DISTORED/DISTORED_max;
end
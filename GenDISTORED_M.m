function DISTORED=GenDISTORED_M(CALC,IRFI_sub,comparerange,TOTSIG)
A=0;
DSUM=0;

DISTORED=zeros(comparerange+1,1);
for ii=1:comparerange+1

    DISTORED(ii,1)=sum(CALC(1:ii,1).*flipud(IRFI_sub(1:ii,1)));
DSUM=DSUM+DISTORED(ii,1);
end

DISTORED=DISTORED.*(TOTSIG/DSUM);

end
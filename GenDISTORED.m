function DISTORED=GenDISTORED(CALC,IRFI_sub,comparerange,TOTSIG)
A=0;
DSUM=0;

%DISTORED=zeros(comparerange+1,1)
for ii=1:comparerange+1

for j=1:ii

    A=A+CALC(j,1)*IRFI_sub(ii-j+1,1);
      
end
DISTORED(ii,1)=A;
DSUM=DSUM+A;
A=0;
end

DISTORED=DISTORED.*(TOTSIG/DSUM);

end
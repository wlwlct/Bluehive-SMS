%Base on the assumpution and my understanding, it is very important to get
%the short wavelength and long wavelength separate, and detect the
%lifetime. Because I think all the prediction we made is based on long
%chain, blue emitter has shorter lifetime, otherwise...how to you explain
%destroy chromophore one by one???

%% This part could be used to generate wavelength change with exp going on
function [ccdt_wavelength,ccdt]=Genccdwl(A)

x_start=A(1,1);
x_end=A(100,1);
%wavelength=A(1:100,1);
for i=1:100
wavelength(i,1)=i;
end

% PART II: Reorganize the CCD data
[nx,ny] = size(A);  
nt = round(nx/100);  
colsum = zeros(1,26);   
% add up the columns to see which is largest
for j = 2:26;
    for i = 101:200;     % because data with light on begin at the second
                        % CCD opening (rows 1-100 are the dark background)
    colsum(1,j) = colsum(1,j) + A(i,j);
    end
end
[maxcol,jmax]=max(colsum);

% array ccdt will have 100 rows of spectral data at nt times
% where each time is a column

ccdt = zeros(100, nt);

% populate the CCDT matrix from the A matrix

for j = 1:nt;
    for i = 1:100;
        k = i + (j-1)*100;   
        ccdt(i,j) = A(k,jmax);
    end
end


    % PART IIC: Dark background subtraction
    
% use the first CCD image for background subtraction

bkg = 0;

for i = 1:100;
    bkg = bkg + ccdt(i,1);
end

bkgavg = bkg/100.;

ccdt = ccdt - bkgavg;    % subtract the average dark pixel from every point



% PART IID: Remove data glitches (despike ccdt(100,nt)
% Insist that values in the CCD data matrix ccdt do not have
% neighboring pixels with more than 3x the intensity and values
% above the average across the third spectrum (choice of third arbitrary) 

av3 = 0;         % average of the third spectrum 

for i = 1:100;
    av3 = av3 + ccdt(i,3);
end
av3 = av3/100.;

for j = 2:nt;   % starts at second time since first spectrum is background
    for i = 1:99;
        if (ccdt(i,j) > 3*ccdt(i+1,j)) && (ccdt(i,j) > av3);
            ccdt(i,j) = 0;  % set spike equal to zero
        end
    end
end

[zong,heng]=size(ccdt);

%mesh(ccdt(20:end,1:end));

for i=2:heng
 
ccdt_wavelength(i,1)=sum(ccdt(20:zong,i).*wavelength(20:end))/sum(ccdt(20:zong,i));

end

end
%%



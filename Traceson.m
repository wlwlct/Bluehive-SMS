%This pfogram works as function,same function as trace, which seperate the
%states. add up different states they calculate.

%This program need to differentiate different energy level and get
%corresponding spectrum and etc....
%This program could assume state has extreme close relation with intensity
%how many state you want based on your model,(or you can automatically generate optimise result )
function [eff,eff_fit,MDL,numst,current_state]=Traceson(countsrate,codefolder)
%%
%clear A adjust breaks columns countsrate current_state delta Del eff eff_fit elimate filename it_endpoint fit_leng fit_startpoint G groups i Ij l leng MDL MDL_min n n_mdl nleng nleng_column numst pathname points pp q rows save_eff_fit sd T1 T2 time Tj uni y  
% PART I: Read in the data


cd(codefolder)
eff=[];
breaks=[];
groups=[];
T1=numel(eff)+1;

%
%This is the part try to use w1_noise detect...
%This part, by using sigma, we could get the medium value of the result
%%Estimate the standard deviation of noise from wavelet at scale =1
delta=diff(countsrate);
y=abs(delta);
y=sort(y);
sd=y(round(0.682*numel(delta)));
%estimate the noise level
sd=sd/1.4;%median absolute deviation

%
%Change_point_wavelet 需要两个variable
%Combine the idea of Haar wavelet and change point method
points=change_point_detection(countsrate);
eff=[eff,countsrate];
T2=numel(eff);
breaks(end+1)=T2;
groups=cat(2,groups,[T1,points+T1;points+T1-1,T2]);
sd=max(sd);

%

%step2 and 3: clustering the segments and calculate MDL
%clusterig_GCP
[G, Ij, Tj] = clustering_GCP(eff, groups);
G = G(end:-1:1);% flip the G
n_mdl = min(30, numel(G));% calculate up to 30 states
MDL = zeros(1,n_mdl);
eff_fit = zeros(n_mdl, numel(eff));
for i = 1:n_mdl
    [MDL(i), eff_fit(i,:)] = MDL_piecewise(Ij, Tj, G(i), eff, groups, sd, breaks);
end


%Following program will be modified from original one.
save_eff_fit=eff_fit;%saving original data...

%
%There was some spikes, elimit spike based on estimate how long it might exist

%Need current state first
numst=1;
MDL_min=min(MDL);
Del=max(MDL)-min(MDL);
[~,pp]=find(MDL<=(MDL_min+0.1*Del));
q=min(pp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%This
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changepoints

    current_state = max(2,q);
    numst=current_state;
    
%In matrix called nleng, we will write down the number in the first row and
%how many points exist in the second row
nleng=zeros(2,1);
[rows,columns]=size(eff_fit);
i=1;n=1;l=1;leng=0;elimate=3;

%Before any calculation, modify spike first. If number less than 10 points I 
%would make them equal to the former number. because I need at least 1 second
%to integrate, and the exp setting is always 0.1s.
 while (i)<=columns
    while eff_fit(numst,i)==eff_fit(numst,i+1)
    i=i+1;
    leng=leng+1;
    if i>=columns
       break 
    end
    end
    nleng(1,n)=eff_fit(numst,i);
    nleng(2,l)=leng+1;
    leng=0;
    i=i+1; n=n+1;l=l+1;
    end

    %Finding the length shorter than certain range, then make the fit intensity
    %equal to the last fit intensity.
    
    [~,nleng_column] = find(nleng(2,:)<=elimate);
    for adjust=1:length(nleng_column)
if nleng_column~=1
    fit_endpoint = sum(nleng(2,1:nleng_column(1,adjust)-1));
    fit_startpoint = sum(nleng(2,1:nleng_column(1,adjust)-1))+1;
    fit_leng=nleng(2,nleng_column(1,adjust));
    eff_fit(numst,(fit_startpoint:fit_endpoint+fit_leng))= eff_fit(numst,fit_endpoint);
end   
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%This part for keeping four states in the diagram
%This part is delete from the current program, I would say what we get is
%what we need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   





end


    



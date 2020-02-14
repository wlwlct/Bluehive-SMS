%GenDISTORED is slow, update to GENDISTRED_M on Feb/19/2018
%Instead of using 0, use -1.

%Only responsible for generate histogram in possible region, not possible
%for generate determin startponit.etc.

%"comparerange" need to be large enough to include all instrumental
%response inside.

%The important part of this program is that it should should distinguish
%background and signal(by find change points),on the other hand, both laser
%light itself or fluorescene need to be normalized,otherwise how do you
%compare with each other...(This is done by the over the max value...)

%Test if the shutter is already closed, and you only get signel from 15
function [how,delta,fitting]=PTU_lifetimefitson_CHM232(Fluo_dtime,tau_min,tau_max,tau_inter,IRFI_hisdtime,IRFI_resolution)
if length(Fluo_dtime)<1000
    how=-1;delta=-1;fitting=-1;return
end

DISTORED=[];
detect_range=[];
%Part I: import instrument response function,called RESP and normalize it
Fluo_max=max(Fluo_dtime);
IRF_data_length=length(IRFI_hisdtime);
IRFI_bg=IRFI_hisdtime(end-1100:end-100);
% %%%Check point%%%%
% figure
% plot(IRFI_bg)
% %%%Check point%%%%
IRFI_bg=sum(IRFI_hisdtime(end-1100:end-100))/1000;
comparerange_CALC=tau_max*10/IRFI_resolution;
comparerange_CALC=max([comparerange_CALC,Fluo_max,IRF_data_length]);
cpst_n=comparerange_CALC-IRF_data_length;
cpst=ones(cpst_n,1)*IRFI_bg;
IRFI_hisdtime=cat(1,IRFI_hisdtime,cpst);%No need to histogram, only need to plot
disp('Get the IRFI back ground')

%Part II: generate CALC:exponential decay in certain range 
tau=400/IRFI_resolution;
[CALC,TOTSIG]=GenCALC_M_PTU(tau,comparerange_CALC);

%Part III: Use above knoledge to generate possible lineshape of
%concolution,called DISTORED. The following part is calculated data.
DISTORED=GenDISTORED_M_PTU(CALC,IRFI_hisdtime,TOTSIG,comparerange_CALC);

%%%Check point%%%%
% figure
% plot(CALC)
% hold on 
% plot(DISTORED)
% hold off
%%%Check point%%%%

%Part IV: This part use the fluorescence to matching the calculated data. I
%guess the main difficulty is the offset. The deviation could be calculated
%For recording DISTORED and Fluo_sub
fitting=[];
%%%
if length(Fluo_dtime)==0
   how=-1;delta=-1;fitting=-1;return
end
[Fluo_N,Fluo_edge]=histcounts(Fluo_dtime,1:1:max(Fluo_dtime));
%Fluo_N_fit=sgolayfilt(Fluo_N,3,111);
Fluo_N_fit=smoothdata(Fluo_N,'sgolay',40);

    %%%Check point%%%%
%     figure
%     plot(Fluo_N_fit);
%     hold on 
%     plot(Fluo_N);    
    %%%Check Point%%%%

Fluo_N_fit_bg=sum(Fluo_N_fit(1,end-1100:end-100))/1000;
Fluo_N=Fluo_N-Fluo_N_fit_bg;

% %%%Check point%%%%
% figure
% histogram(Fluo_dtime,Fluo_edg)
% hold on
% plot(Fluo_N)
% %%%Check point%%%%

%Part IV.1: The Fluo_N has limited data which is useless for cauculating
if length(find(Fluo_N>=5))==0
    how=-1;delta=-1;fitting=-1;
   return
end

if length(find(Fluo_N~=0))<=700
     how=-1;delta=-1;fitting=-1;
   return
end

Fluo_N_fit=smoothdata(Fluo_N,'sgolay',40);
[Fluo_N_fit_max,Fluo_N_fit_max_position]=max(Fluo_N_fit);

[DISTORED_max,DISTORED_max_position]=max(DISTORED);
position_difference=Fluo_N_fit_max_position-DISTORED_max_position;
if position_difference>0
   Fluo_sub=cat(2,Fluo_N(1,position_difference:end),Fluo_N(1,1:position_difference-1));
else
    Fluo_sub=cat(2,Fluo_N(1,end+position_difference+1:end),Fluo_N(1,1:end+position_difference));
end
Fluo_sub_length=length(Fluo_sub);
if length(DISTORED)>=Fluo_sub_length
   DISTORED=DISTORED(1:Fluo_sub_length,1); 
end
Fluo_sub=Fluo_sub./Fluo_N_fit_max;

    %%%Check point%%%%
%     figure
%     plot(transpose(Fluo_sub));
%     hold on 
%     plot(DISTORED)
%     hold off
    
    %%%Check Point%%%%

%Part V: Find the parameter to determine the how good is the fitting.
ddd=0;
lleenngg=length(tau_min/IRFI_resolution:tau_inter/IRFI_resolution:tau_max/IRFI_resolution);
what.residue=zeros(length(Fluo_sub),lleenngg);
what.S=zeros(lleenngg,2);
for tau=tau_min/IRFI_resolution:tau_inter/IRFI_resolution:tau_max/IRFI_resolution  
    [CALC,TOTSIG]=GenCALC_M_PTU(tau,comparerange_CALC);
    DISTORED=GenDISTORED_M_PTU(CALC,IRFI_hisdtime,TOTSIG,comparerange_CALC);   
    [DISTORED_max,DISTORED_max_position]=max(DISTORED);
    position_difference=Fluo_N_fit_max_position-DISTORED_max_position;
    if position_difference>0
        Fluo_sub=cat(2,Fluo_N(1,position_difference:end),Fluo_N(1,1:position_difference-1));
    else
        Fluo_sub=cat(2,Fluo_N(1,end+position_difference+1:end),Fluo_N(1,1:end+position_difference));
    end
    Fluo_sub=Fluo_sub./Fluo_N_fit_max;  
        
    check_n=1;
    small_chisqr=1;
    Fluo_sub_copy=Fluo_sub;
    %shift Fluo left and right to get better fitting position
    while check_n<50
        check_n=check_n+1;
        sqrr=DISTORED(1:DISTORED_max_position+200,1)-transpose(Fluo_sub_copy(1,1:DISTORED_max_position+200));
        CHIsum=sum(sqrr.^2./DISTORED(1:DISTORED_max_position+200,1));
        chisqr=chi2cdf(CHIsum,DISTORED_max_position+200);
        if small_chisqr>chisqr
            small_chisqr=chisqr;
            Fluo_sub=Fluo_sub_copy;
        end
        if chisqr<10^(-11)
            break
        else
            shift=2;
            maxlength=40*shift;
            if check_n==2
                Fluo_sub_copy=cat(2,Fluo_sub_copy(1,maxlength:end),Fluo_sub_copy(1,1:maxlength-1));
            else
                Fluo_sub_copy=cat(2,Fluo_sub_copy(1,end-shift+1:end),Fluo_sub_copy(1,1:end-shift));
            end 
        end
   end
    
    Fluo_sub_length=length(Fluo_sub);
    if length(DISTORED)>=Fluo_sub_length
        DISTORED=DISTORED(1:Fluo_sub_length,1); 
    end

    ddd=ddd+1;
    what.S(ddd,1)=sum((transpose(Fluo_sub)-DISTORED).^2)/Fluo_sub_length;%This is S value
    what.S(ddd,2)=tau*IRFI_resolution;%Unit of tau is bin number and Unit of other tau is ps.
    what.ft(ddd).tau=cat(2,transpose(Fluo_sub),DISTORED);
    
    %%%Check point%%%%
%     figure
%     plot(transpose(Fluo_sub));
%     hold on 
%     plot(DISTORED)
%     hold off
%     
    %%%Check Point%%%%
    what.residue(:,ddd)=(DISTORED-transpose(Fluo_sub));   
end

if size(what)~=0
    [~,index]=mink(what.S(:,1),ceil(length(what.S(:,1))*0.1));%only seleect 10% of the fitting
    for dddd=1:length(index)
        fitting.S(dddd,1)=what.S(index(dddd,1),1);%record S
        fitting.S(dddd,2)=what.S(index(dddd,1),2); %record tau
        fitting.ft(dddd).tau=what.ft(index(dddd)).tau;%contain distored and fluo_sub
    end
    fitting.residue=what.residue(:,index(:,1));
    delta=max(fitting.S(:,2))-min(fitting.S(:,2));
    how(:,1)=fitting.S(:,1);%Record S in first column, and tau in second column.
    how(:,2)=fitting.S(:,2);
else
    how=-1;fitting=-1;delta=-1;
end
end

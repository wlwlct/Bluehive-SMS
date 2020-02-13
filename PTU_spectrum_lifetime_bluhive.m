%02/12/2020
clear all
clc;
%%This part for generate how wavelength change with experiment time, assume
%%integration time is 1 sec.
%%
%PartI: user define
%User define. This part is for chose the direction of apd file, then automatically run files
codefolder=pwd;
solvent='F8T2 Chloroform 2kda ';
%date='$date'%basically the folder name
date='02022020'
%datafolder=['../' date '.$foldernumber'];
datafolder=['C:\Users\Livi\Documents\Results\' date];
IRFI_location='C:\Users\Livi\Documents\Results\02022020\IRFI_8ps_fakegaussian7_channel'
cd([datafolder '/apd']);
%deal with APD file
allnames=struct2cell(dir('*.mat'));
len=length(allnames(1,:));
inttime=1;%This is time in sec
deadtime=0.079;%This is time in sec
timetrace_resolution=0.1*10^9;%The resolution for setting the bin size in time trace in ns
tau_max=2500; % in unit of picoescond
tau_min=30;%in unit of pico second
tau_inter=10;%between lowest and highest tau, you won't check one by one....
plot_graph='no';
dataset_Y='yes';
channel_box=[0];%must put in a row
wavelengthstart=1;
disp('Finish User define\n')

% Part II. start the calculation
for len_i=1:len;
    clearvars -except len_i date k len allnames inttime deadtime plot_graph dataset_Y timetrace_resolution codefolder channel_box IRFI_location tau_max tau_min tau_inter datafolder solvent wavelengthstart
    try
        clear name
        name=char(allnames(1,len_i));
        name=name(1:length(name)-4);
        cd(codefolder);
        %import CCD
        try
            ccdt=importdata([datafolder '/ccd/' 'ccdt' name '.mat']);
            expwl=importdata([datafolder '/ccd/' 'ccdt_wavelength' name '.mat']);
        catch
            fprintf('No CCD file avaliable for file: %s\n',name)
        end
        %Import APD file.
        try
            apd_file=[datafolder '/apd/' name '.mat'];
        catch
            fprintf('No APD file avaliable for file: %s\n. Not suppose to happen\n',name)
        end
        
        try
            for channel_choice=1:numel(channel_box)
                cd(codefolder);
                channel=channel_box(channel_choice,1);
                %Part II.1 Time trace of APD
                [apddata,apddataresolution]=PTUim(apd_file);
                disp('Finish import apdfile\n')
                datasource=GetDandABS(apddata,channel,'M');
                absolutetime=datasource(:,1);
                dtime=datasource(:,2);
                %Deserve a for loop for different data sets
                timetrace=Gen_timetrace(absolutetime,timetrace_resolution+deadtime*(10^9));
                disp('Finish generate time trace')
                time=timetrace(:,1);
                %PartII.2 Modified curve calculation
                countsrate=timetrace(:,2);%Counts in 10ms
                [eff,eff_fit,MDL,numst,current_state]=Traceson(countsrate,codefolder);%this for seperate into different states
                disp('Finish generate modified time trace')
                %Generate the place the change occur,the last element of former segment, not the last segment;
                stage_changepts=transpose(find(diff(eff_fit(numst,:))~=0));
                stage_start=time(stage_changepts,1);%Generate time corresponding to stage change              
                stage_start=[0;stage_start;time(end)];%make the unit from sec to nanosec,and start with time zero
                disp('Finish Generate change point with absolute time')
                %This is for dissect data to 1 sec time range, or several time range when points are not enough.
                ttstart=absolutetime(1,1);%The unit is in nanosecond
                ttend=absolutetime(end,1);
                exptime = ttend;
                perstartt = ttstart;
                
                perfect_prep=[0:1:ceil((ttend-ttstart)/((inttime+deadtime)*10^9))];
                perfecttime=perstartt+perfect_prep.*((inttime+deadtime)*10^9);%matrix for the desired timeline
                perfecttime=perfecttime(perfecttime<exptime-inttime*10^9);           
                disp('Finish Generate perfect time line')
                %%This part to seperate the perfectime in different state, so we can use
                %seperate different row range in different state.Based on the state,
                %if unique 2 states, then calculate fist 1/2, and last 1/3,if more than
                %unique 2 states, then just calculate based on each second first, if you
                %couldn't get any lifetime based on one second, then just add up following
                %rowrange, untill we can get a proper lifetime. If we couln't
                %(unlikely),just let it go.

                %
                rowrange = [];
                stagestart_n=2;
                strow=0;
                lengp=length(perfecttime);

                for n = 1:lengp
                    [~,srn]=min(abs(absolutetime-perfecttime(1,n)));
                    [~,ern]=min(abs(absolutetime-perfecttime(1,n)-inttime*10^9));
                    
                    %This part use two 'if' to put data into different state...
                    if stagestart_n<=length(stage_start)
                        if stage_start(stagestart_n-1,1)<=perfecttime(1,n) && perfecttime(1,n)<=stage_start(stagestart_n,1)
                            strow = strow+1;
                            rowrange(stagestart_n-1).rr([1,2],strow)=[srn;ern];
                        elseif perfecttime(1,n)>stage_start(stagestart_n,1)
                            strow=1;
                            stagestart_n=stagestart_n+1;
                            rowrange(stagestart_n-1).rr([1,2],strow)=[srn;ern];
                        end        
                    else
                        strow = strow+1;
                        rowrange(stagestart_n-1).rr([1,2],strow)=[srn;ern];
                    end
                end
                
                %%%
                for n = 1:lengp
                    [~,perfect_rowrange(1,n)]=min(abs(absolutetime-perfecttime(1,n)));
                    [~,perfect_rowrange(2,n)]=min(abs(absolutetime-perfecttime(1,n)-inttime*10^9));
                end
                %%%
                
if length(rowrange(1).rr)==0;
    rowrange(1)=[];
end
disp('Finish use perfet time line to generate row range')
cd(codefolder)
%
%This is the part for generate lifetime for each part
IRFI_sub_file=strcat([IRFI_location num2str(channel) '.mat']);
[IRFI_hisdtime,IRFI_resolution]=PTUim(IRFI_sub_file);
disp('Finish load IRFI file')
if IRFI_resolution~=apddataresolution
disp('No need to calculate lifetime, IRF and Fluorescence have different resolution')
end
    
%Part 0: User Define
%Part0.2 program required
lf=[];
%For generating picture purpose, I will introduce variable 'js' to check
%what is the second number.
%%%for name
js=0;
%%%
for n=1:length(rowrange)
    [~,rowrange_heng]=size(rowrange(n).rr);
    zeroplace(n).zp(1,1)=0;
    ii=0;
for i=1:rowrange_heng
    %%%For name
    js=js+1;
    %%%
    Intensity(js)=length(rowrange(n).rr(1,i):rowrange(n).rr(2,i));
disp('Finish Intensity part');
    Fluo_dtime=dtime(rowrange(n).rr(1,i):rowrange(n).rr(2,i),1);
    disp('Finish dtime part')
    
    
[Live(n).lifetime(i).ll,~,fitting(n).fit(i).ft]=PTU_lifetimefitson_CHM232(Fluo_dtime,tau_min,tau_max,tau_inter,IRFI_hisdtime,IRFI_resolution);%By using this, need to find your own range....

if Live(n).lifetime(i).ll == 0
    %find the place equal to zero in each state.
        ii=ii+1;
    zeroplace(n).zp(ii,1)=i;
end
end
fprintf('Finish generate lifetime related information in stage %d',n);
%
%Also you need to check which part is continuous zero. I check by minus
%former value. I also put different set of continuous part in 'matrix
%conti',column direction is different set



%If zero place is inside, one of the problem is that zero won't produce
%unless the next value shows up, I think...
if length(zeroplace(n).zp)>=2;j=1;k=1;q=1;
for i=1:length(zeroplace(n).zp)-1
        if zeroplace(n).zp(i,1)==zeroplace(n).zp(i+1,1)-1
        q=0;
        conti(n).co(j,k)=zeroplace(n).zp(i,1);
        k=k+1;
            if i==length(zeroplace(n).zp)-1
                conti(n).co(j,k)=zeroplace(n).zp(i+1,1);
            end
        
        else
            if q~=1
            conti(n).co(j,k)=zeroplace(n).zp(i,1);
            j=j+1;% j for each continuous 1 put in one row,change different row for different set of       
            q=1;k=1;
                if i==length(zeroplace(n).zp)-1
                conti(n).co(j,k)=zeroplace(n).zp(i+1,1);
                end
            end
        end
        
end        
end
end
disp('Finish check place lifetime equal 0')
%%
%This part is another for loop for calculate rearranged data.
%If the zero is continuous, add up nearby two set of data, recalculate lifetime again;if the zero is not
    %continuous,give up there might be some extra problem that we need to
    %solve. Similar to Lifetime calculation

    
    %Avoid the confucion and easier the cauculation, if the 0 more than 
    %js_n=0;
for n=1:length(rowrange)
    [conti_zong,~]=size(conti(n).co);
  
    Fluo_dtime=[];
    if numel(conti(n).co)~=0
    for iii=1:conti_zong
      zerolength=length(find(conti(n).co(iii,:)>0));
                if zerolength<30;
                    Fluo_dtime=dtime(rowrange(n).rr(1,conti(n).co(iii,1)):rowrange(n).rr(2,conti(n).co(iii,1)),1);
                    for k=1:zerolength-1;
                    Fluo_dtime=cat(1,Fluo_dtime,dtime(rowrange(n).rr(1,conti(n).co(iii,1+k)):rowrange(n).rr(2,conti(n).co(iii,1+k)),1));                  
                    end
                    %%%For name
                    %js=js_n+conti(n).co(iii,1);               
                   %%%
                    [Live(n).lifetime(conti(n).co(iii,1)).ll, ~,fitting(n).fit(conti(n).co(iii,1)).ft]=PTU_lifetimefitson_CHM232(Fluo_dtime,tau_min,tau_max,tau_inter,IRFI_hisdtime,IRFI_resolution);                    for k=1:zerolength-1
                    Live(n).lifetime(conti(n).co(iii,1+k)).ll=Live(n).lifetime(conti(n).co(iii,1)).ll;
                    fitting(n).fit(conti(n).co(iii,1+k)).ft=fitting(n).fit(conti(n).co(iii,1)).ft;
                    end
                    
                elseif zerolength>=30
      %Then we need to calculate from both side, then compare with each
      %other, basically we need to seperate the zeros into 3 parts, with
      %intersection of 10s, then we will calculate from both side, then
      %lets see what's the proble lifetime for both side.
                    callength=floor((zerolength-10)/2)-1;%with minus one to eliminate the chance of over exceed the length
                    %for claculate and write in the first half
                    Fluo_dtime=dtime(rowrange(n).rr(1,conti(n).co(iii,1)):rowrange(n).rr(2,conti(n).co(iii,1)),1);
                    for k=1:callength-1
                    Fluo_dtime=cat(1,Fluo_dtime,dtime(rowrange(n).rr(1,conti(n).co(iii,1+k)):rowrange(n).rr(2,conti(n).co(iii,1+k)),1));                  
                    end
                     %%%For name
                     %js=js_n+conti(n).co(iii,1); 
                     %%%
                    [Live(n).lifetime(conti(n).co(iii,1)).ll, ~,fitting(n).fit(conti(n).co(iii,1)).ft]=PTU_lifetimefitson_CHM232(Fluo_dtime,tau_min,tau_max,tau_inter,IRFI_hisdtime,IRFI_resolution);
                     for k=1:callength-1
                        Live(n).lifetime(conti(n).co(iii,1+k)).ll=Live(n).lifetime(conti(n).co(iii,1)).ll;
                        fitting(n).fit(conti(n).co(iii,1+k)).ft=fitting(n).fit(conti(n).co(iii,1)).ft;
                    end
                    
                    %for calculate and write in the second half 
                    k=callength+10;
                    kk=callength+10;
                    Fluo_dtime=dtime(rowrange(n).rr(1,conti(n).co(iii,k)):rowrange(n).rr(2,conti(n).co(iii,k)),1);
                    for k=callength+10:zerolength-1
                    Fluo_dtime=cat(1,Fluo_dtime,dtime(rowrange(n).rr(1,conti(n).co(iii,1+k)):rowrange(n).rr(2,conti(n).co(iii,1+k)),1));
                    end
                     %%%For name
                     %js=js_n+conti(n).co(iii,k);
                     %%%
                    [Live(n).lifetime(conti(n).co(iii,kk)).ll, ~,fitting(n).fit(conti(n).co(iii,kk)).ft]=PTU_lifetimefitson_CHM232(Fluo_dtime,tau_min,tau_max,tau_inter,IRFI_hisdtime,IRFI_resolution);
                    
                    for k=callength+10:zerolength-1
                        Live(n).lifetime(conti(n).co(iii,k+1)).ll=Live(n).lifetime(conti(n).co(iii,kk)).ll;
                        fitting(n).fit(conti(n).co(iii,k+1)).ft=fitting(n).fit(conti(n).co(iii,kk)).ft;
                    end
                end      
        %also later would add on how to deal with different stage, two stage or
    %two more, we will treat them differently.
    end
    end   
    
    %%%For name
     %   [~,rowrange_heng]=size(rowrange(n).rr);
    %js_n=js_n+rowrange_heng;
    %%%
    
end
disp('Finish recombination of 0 places and recalculate')

%%
% this part could transfer the ensemble data from different segments to one
% matrix with lifetime for each second.
i=0;
for n=1:length(Live)
    [~,lifetime_heng]=size(Live(n).lifetime);
    for ii=1:lifetime_heng
        [lf_min,lf_min_position]=min(Live(n).lifetime(ii).ll(:,1));
    if lf_min~=0
    lf(i+ii,1)=lf_min;%for If the first column is the min S value,second is corresponsing lifetime, third is minimum lifetime, forth corresponding S value; fifth is maximum lifetime; 6th is correponding S value.
    lf(i+ii,2)=Live(n).lifetime(ii).ll(lf_min_position,2);%lifetime for mininum S value
    [lf(i+ii,5),min_lifetime_position]=min(Live(n).lifetime(ii).ll(:,2));%minimum lifetime 
    lf(i+ii,3)=max(Live(n).lifetime(ii).ll(min_lifetime_position,1));%S value corresponding to minimum lifetime
    [lf(i+ii,6),max_lifetime_position]=max(Live(n).lifetime(ii).ll(:,2));%maximum lifetime
    lf(i+ii,4)=max(Live(n).lifetime(ii).ll(max_lifetime_position,1));%S value corresponding to maximum lifetime
    end
    end
    i=lifetime_heng+i;
end

disp('Finish find S related value');
%%
%plot the graph
if length(lf)~=0

[lf_v,~]=size(lf);
[wl_v,~]=size(expwl);
%plot
    if lf_v>=wl_v
    %lf_max=max(lf(1:wl_v-1,2));%max value is for normalization, may bring some
    %erro, not use here for now
    %expwl_max=max(expwl(2:end,1));
    lfneg=lf(1:wl_v-1,5)-lf(1:wl_v-1,2);
    lfpos=lf(1:wl_v-1,6)-lf(1:wl_v-1,2);
    lfS=lf(1:wl_v-1,1);
    lfwl=cat(2,lf(1:(wl_v-1),2),expwl(2:end,1));%minus one is due to first element of expwl is 0. in order to make the dimension match each other.
    ccdintensity=sum(ccdt(:,3:end));
disp('Finish if to cut data \n')
else
    %lf_max=max(lf(:,2));
    %expwl_max=max(expwl(2:lf_v,1));
    lf=cat(1,lf,zeros(wl_v-lf_v-1,6));
    lfneg=lf(:,5)-lf(:,2);
    lfpos=lf(:,6)-lf(:,2);
    lfS=lf(:,1);
    lfwl=cat(2,lf(:,2),expwl(2:end,1));% same as minus one
    ccdintensity=sum(ccdt(:,3:end));
disp('Finish if to cut data \n')
    end
    
%check combineplace, write new ccd and average ccd part 
    xupper=length(expwl(2:end,1));
[allsecs,newconti]=con2sec(rowrange,conti,xupper);
newccdt=ccdt;
newexpwl=expwl;

newconti_leng=length(newconti);
for newconti_i=1:1:newconti_leng
newcontico_leng=length(newconti(newconti_i).co);
for newcontico_i=1:1:newcontico_leng
    
    if ~isempty(newccdt(:,newconti(newconti_i).co(newcontico_i).subco))
newccdt_combination=sum(newccdt(:,newconti(newconti_i).co(newcontico_i).subco(1,1)+2:...
    newconti(newconti_i).co(newcontico_i).subco(1,end)+2),2);
avewavelength=sum(newccdt_combination(wavelengthstart:end,1).*newccdt(wavelengthstart:end,1))/sum(newccdt_combination(wavelengthstart:end,1));

subco_leng=length(newconti(newconti_i).co(newcontico_i).subco(1,:));
for subco_i=1:1:subco_leng

    newccdt(:,newconti(newconti_i).co(newcontico_i).subco(1,subco_i)+2)=newccdt_combination;
    newexpwl(newconti(newconti_i).co(newcontico_i).subco(1,subco_i)+1,1)=avewavelength;
    
end


    end



end

end



disp('Successful with newconti')
  %This is for generate the useful dataset
  if strcmp(dataset_Y,'yes')==1
     dataset.time_trace=cat(2,timetrace(:,1),timetrace(:,2).*(0.01*10^9)/timetrace_resolution);
     dataset.side=struct('eff',eff*((0.01*10^9)/timetrace_resolution),'eff_fit',eff_fit*((0.01*10^9)/timetrace_resolution),'MDL',MDL,'numst',numst,'current_state',current_state,'ABStime_x',time-min(time));
     dataset.ccdt=ccdt(:,1:end);
     dataset.newccdt=newccdt(:,1:end);
     dataset.scatterplot.lifetime=cat(2,lfS,lfwl(:,1),lfneg,lfpos);
     dataset.scatterplot.spectrum=lfwl(:,2);
     dataset.scatterplot.newspectrum=newexpwl(2:end,1);
     dataset.scatterplot.intensity=cat(1,Intensity(1,1:length(lfwl(:,1))),ccdintensity);
     dataset.fitting=fitting;
     dataset.rowrange=rowrange;
     dataset.newconti=newconti;
     dataset.allsecs=allsecs;
      
        try
    cd ([datafolder '/dataset intermediates' '/' num2str(channel)])
        catch
    cd(datafolder)
    %mkdir dataset intermediates
    mkdir([datafolder '/dataset intermediates' '/' num2str(channel)])
    cd ([datafolder '/dataset intermediates' '/' num2str(channel)])
        end
    save([solvent date 'dataset' name '.mat'],'dataset')                
  end
    %
   disp('Finish save dataset')


if strcmp(plot_graph,'yes')==1
    
figure
hold off
subplot(3,4,[7 8 11 12]);
title('lifetime and spectrum')
x=1:length(lfwl(:,1));
%scatter3(x,lfwl(:,1),lf(x,1),'o');
errorbar(x,lfwl(:,1),lfneg,lfpos, 'o');
text(x,lfwl(:,1),lf(x,1)-1,'\ast','HorizontalAlignment','center')
zlim([-0.15 0.15]);


hold on 
plot(lfwl(:,2),'o');
legend('lifetime','wavelength (pixel*8)')%using normalized picture will miss some important points, such as in raster scan, the dark position seems to be blue emission.
end
else
    if strcmp(plot_graph,'yes')==1
    subplot(3,4,[7 8 11 12]);
    title('spectrum')
    plot(expwl.*8,'o');
    xlabel('Lifetime not work for this file...sad...')
    %print=['Not working for this file....sad....']
    end
end
%
%This part plot time trace in the subplot of (3,4,[1 2])(3,4,3)(3,4,4)
if strcmp(plot_graph,'yes')
    numst = min(numst, numel(MDL));
    current_state = numst;
    subplot(3,4,3);
    %state plot is for histogram
    cd(codefolder);
    state_plot(numst,eff,eff_fit,MDL);
    title('Intensity horizontal histogram')
    
    subplot(3,4,[1 2]);
    slider_plot(eff,eff_fit,current_state);
    title(strcat('Time Trace ', ' Choosing ',num2str(numst),' states'))
    %rowrange
    subplot(3,4,4);
    spider_plot(MDL,numst);
    title('MDL vs. States')

%%
%This part for generate ccd spectrum
subplot(3,4,[5 6 9 10]);
hold
mesh(ccdt);title('spectrum');
%%
%This is the part for saving graph
try
cd ([datafolder '/Figure intermediates' '/' num2str(channel)])
catch
    cd([datafolder])
    mkdir Figure intermediates
    cd ([datafolder '/Figure intermediates' '/' num2str(channel)])
end
savefig([solvent date name '.fig'])
close all
end
end
catch
fprintf('file: %s channel %d fail',name,channel_choice)
end


%trycatch for len loop
catch
    fprintf('file: %s fail',name)
end
end

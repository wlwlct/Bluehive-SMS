function Fluo_sub=GenFluo_sub(Fluo_N,Fluo_startpoint,Fluo_max,comparerange,DISTORED)


Fluo_sub=transpose(Fluo_N(1,Fluo_startpoint:(Fluo_startpoint+comparerange)));
%%%%test
Fluo_subb=sgolayfilt(Fluo_sub,3,21);%This for peak choosing compensate
Fluo_max=max(Fluo_subb);

%%%%%%%%%%%%%%%%%%%%%%%%%test for compensate the fluffy max value may lead
%%%%%%%%%%%%%%%%%%%%%%%%%to some error


Fluo_sub=Fluo_sub./Fluo_max;
%plot(Fluo_sub)
%hold on
%%%test
%Fluo_sub=GenFluo_sub(Fluo_sub,Fluo_max,comparerange,DISTORED);


%make the fluorescence graph maximum matching calculated maximum position
%and same height ratio

%Find the height and location of calculated graph
[DISTORED_max,DISTORED_max_row]=max(transpose(DISTORED));
%Change the size of fluorescnce data %change the range to ratio



Fluo_sub=Fluo_sub.*DISTORED_max;

%change the off set of the fluorescence.
[~,Fluo_max_row]=max(Fluo_sub);
Fluo_new_startpoint=Fluo_startpoint+(Fluo_max_row-DISTORED_max_row);
Fluo_sub=transpose(Fluo_N(1,Fluo_new_startpoint:(Fluo_new_startpoint+comparerange)));
Fluo_sub=Fluo_sub./Fluo_max.*DISTORED_max;



end
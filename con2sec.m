  function [allseconds,new]=con2sec(rowrange,conti,xupper)
length_conti=length(conti);
new=struct('co',cell(length_conti,1));
allseconds=[];
for i=1:1:length_conti
        new(i).co=[];
        if ~isempty(conti(i).co)
            
            
    co_height=length(conti(i).co(:,1));
    new(i).co=struct('subco',cell(co_height,1))
    dd=1;
    for ii=1:1:co_height
        buffer=conti(i).co(ii,(conti(i).co(ii,:)~=0));
       iii=i;
       while i>=2 && iii>=2
           buffer=buffer+length(rowrange(iii-1).rr(1,:));
           iii=iii-1;
       end
       buffer(buffer>xupper)=[];    
       
       if length(buffer)>=30
           callength=floor((length(buffer)-10)/2)-1;
           subbuffer{1}=buffer(1,1:callength);
           subbuffer{2}=buffer(1,callength+10:end);
       new(i).co(dd).subco=subbuffer{1};
       dd=dd+1;           
       allseconds=cat(2,allseconds,subbuffer{1});
       new(i).co(dd).subco=subbuffer{2};
       dd=dd+1;           
       allseconds=cat(2,allseconds,subbuffer{2});
       else
       new(i).co(dd).subco=buffer;
       dd=dd+1;
       allseconds=cat(2,allseconds,buffer);
       end 
    end
     
        end
    
end

new_i=1;
while new_i<=length(new)
    if isempty(new(new_i).co)
        new(new_i)=[];
    else
        new_i=new_i+1; 
    end
    
end




  end
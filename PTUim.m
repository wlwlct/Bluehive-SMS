function[data_raw,resolution]=PTUim(fullpath)
A=load(fullpath);
try
data_raw=A.PTU3file.data;
catch 
    data_raw=A.PTU3file;
end
if isempty(data_raw)
   error('Attention! The data import is empty!') 
end
    
resolution=A.Resolution*10^12;
end
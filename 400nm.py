import os
import sys
import math


def useTemplate(templateFile, Values,destination, folder= '~/Template/' ):
    """
    Use a template file 
    to create files
    """
    from string import Template
    filein = open( folder+"/"+templateFile )
    src = Template( filein.read() )
    out = src.safe_substitute(Values)
    outName = open(destination+"/"+templateFile,"w+")
    outName.write(out)
    outName.close()






def dd(command):
	os.system(command)
# Input folder that contains all data
base =  sys.argv[1]
if base[-1]=="/":
	base = base[:-1]


print base 
#base = '/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/'

# Task 1 : 

# gather the filename alphabetically
os.chdir(base + "/apd")
filelist = os.popen("ls -1 | sort").read().split("\n") 
# Filter file names
filterlist = []
i =0
while i < len(filelist):
	if filelist[i].split(".")[-1] == "mat":
		filterlist.append(filelist[i])
	i=i+1

print "Filering names"
# Create many folders
setlength = 5
nfiles = len(filterlist) 
nfolder = int(math.ceil(nfiles/float(setlength)))
print "Will create", nfolder
# creating n folder
os.chdir("../") # will create inside base  
n = 0 
while n < nfolder:
	foldername = base.split("/")[-1] + "." + str(n) 
	print "Creating", foldername
	try :
		os.mkdir(foldername) 
		os.mkdir(foldername+"/apd") 
		os.mkdir(foldername+"/ccd" ) 
	except:
		print "Warning", foldername, "Exists"
	n = n + 1

# copy files
print "Folders Created"

basefolder =  base.split("/")[-1] 
print basefolder 
for i in range(nfiles):
	ithfolder = int(math.ceil(i/setlength)) 
	filemat = "apd/" + filterlist[i] 
	fileccd = "ccd/ccdt" + filterlist[i] 
	fileccdwave =  "ccd/ccdt_wavelength" + filterlist[i]
	os.system("cp %s %s.%s/apd/"%(filemat,basefolder,ithfolder) )
	os.system("cp %s %s.%s/ccd/"%(fileccd,basefolder,ithfolder) )
	os.system("cp %s %s.%s/ccd/"%(fileccdwave,basefolder,ithfolder) ) 

print "Files copied to generated folders"

#create sbatch file
#useTemplate("run.sbatch", params, "/home/lwang74" ,"/home/lwang74/code")  
#useTemplate("PTU_spectrum_lifetime_bluhive.m", {'date':'1009'}, "/home/lwang74" ,"/home/lwang74/code")  
for i in range(nfolder):
	date =  basefolder	
	params = { 'date': date,
		   'base': base,
		   'foldernumber': str(i) }
	os.system("mkdir code%s"%(str(i)))
	os.system("cp /scratch/lwang74/PTU_spectrum_lifetime_bluehive/code4003d3d315/* code%s/"%(str(i)))
	useTemplate("PTU_spectrum_lifetime_bluhive.m", params, "code%s"%(str(i)) ,"/scratch/lwang74/PTU_spectrum_lifetime_bluehive/code4003d3d315")  
	useTemplate("run.sbatch", {'jobname':'%s.%s'%(date,str(i)) }, "code%s"%(str(i)) ,"/scratch/lwang74/PTU_spectrum_lifetime_bluehive/code4003d3d315")  
	os.chdir("code%s/"%(str(i)))
	os.system("sbatch run.sbatch")
	os.chdir("../")


a #!/bin/env python

from os import popen
from os import listdir
from os.path import isfile, join
mypath = '/afs/cern.ch/user/t/taroni/eos/cms/store/user/taroni/ZeroBias/turnOnNtuple/170413_154606/'
#print mypath, listdir(mypath)
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#numero iniziale
i=0
popen("mkdir mp_current") 

#inizio ciclo
while i < len(onlyfiles):

    out_file = open("mp_current/inputList"+str(i)+".txt","w")
    j=0
    while j<1:
        if i*1+j < len(onlyfiles) :
            out_file.write("/store/user/taroni/ZeroBias/turnOnNtuple/170413_154606/"+onlyfiles[i*1+j]+"\n")
        j+=1
    out_file.close()
    
    sh = """#!/bin/tcsh -f
    	  set W_DIR = \"/afs/cern.ch/work/t/taroni/private/EndcapFG/src/TurnOn\"
	  cp $W_DIR/makePairs.cc .     
	  cp $W_DIR/mp_current/inputList"""+str(i)+""".txt .
          cp $W_DIR/baseFuncNad.h . 
          cp $W_DIR/readJSONFile.cc . 
          cp $W_DIR/readJSONFile.h . 
           mkdir makePairsDir
            cd $W_DIR
            eval `scramv1 runtime -csh`
            cd -
            g++ -o makePairs makePairs.cc  `root-config --cflags --libs` -L $ROOFITSYS/lib -lRooFit -lRooFitCore
            ./makePairs myid inputList"""+str(i)+""".txt
            cmsStage -f makePairsDir/elepairs_tree.root /store/user/taroni/ZeroBias/turnOnNtuple/makePairs/elepairs_tree_"""+str(i)+""".root
            exit
    """

    #scrive script
    sh_file = open("mp_current/mp-"+str(i)+".sh","w")
    sh_file.write(sh)
    sh_file.close()

    #sottomette script
    popen("cd mp_current") 
    popen("chmod a+x  mp_current/mp-"+str(i)+".sh" )
    popen("cd  mp_current; bsub -q cmscaf1nd mp-"+str(i)+".sh" )
    
    i+=1


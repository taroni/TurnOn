#!/bin/tcsh -f
    	  set W_DIR = "/afs/cern.ch/work/t/taroni/private/EndcapFG/src/TurnOn"
	  cp $W_DIR/selectPairs.cc .     
	  cp $W_DIR/listOfFiles .
          cp $W_DIR/baseFuncNad.h . 
          cp $W_DIR/readJSONFile.cc . 
          cp $W_DIR/readJSONFile.h . 
          mkdir selectPairsDir
          cd $W_DIR
          eval `scramv1 runtime -csh`
          cd -
          g++ -o selectPairs selectPairs.cc `root-config --cflags --libs` -L $ROOFITSYS/lib -lRooFit -lRooFitCore
          ./selectPairs listOfFiles

          exit
    

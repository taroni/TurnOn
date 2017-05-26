JOBID         = ratioScan
CC            = g++ -g
#LIBS          = `root-config --cflags --libs` -L $(ROOFITSYS)/lib -lRooFit -lRooFitCore -I$(ROOFITSYS)/include 
LIBS          = `root-config --cflags --libs` -L $(ROOTSYS)/lib -lRooFit -lRooFitCore -I$(ROOTSYS)/include


SRCS          = FuncCB.cc fitEfficiency.cc 
SRCSERR       = FuncCB.cc FindValueErr.C 
SELECTSRC     = selectPairs.cc
MAKESRC       = makePairs.cc
RESOLUTIONSRC = resolution.cc

OBJS          = FuncCB.o fitEfficiency.o 

FITERROR      =FindValueErr
FIT           = fitEfficiency
SELECT        = selectPairs
MAKEPAIRS     = makePairs
RESOLUTION    = resolution

all: $(FIT) $(FILE_EXISTS)

$(FIT):   $(SRCS) $(SELECT) 
	echo "Linking $(FIT) ..."
	g++  -o $(FIT) $(SRCS) $(LIBS)


$(FITERROR):   $(SRCSERR) 
	echo "Linking $(FITERROR) ..."
	g++  -o $(FITERROR) $(SRCSERR) $(LIBS)

$(SELECT): $(SELECTSRC) $(MAKEPAIRS) 
	g++  -o $(SELECT) $(SELECTSRC) $(LIBS)
	test -d $(JOBID) || make $(MAKEPAIRS)  make runPairs
	test -d $(JOBID)/selectPairsDir || mkdir $(JOBID)/selectPairsDir

$(MAKEPAIRS): $(MAKESRC) 
	g++  -o $(MAKEPAIRS) $(MAKESRC) $(LIBS)
	test -d $(JOBID) || mkdir $(JOBID)
	test -d $(JOBID)/makePairsDir || mkdir $(JOBID)/makePairsDir

$(RESOLUTION): $(RESOLUTIONSRC) 
	g++  -o $(RESOLUTION) $(RESOLUTIONSRC) $(LIBS)
	test -d $(JOBID) || mkdir $(JOBID)
	test -d $(JOBID)/resolutionDir || mkdir $(JOBID)/resolutionDir


runPairs: $(MAKEPAIRS)
	./$(MAKEPAIRS)  $(JOBID)
runSelect: $(SELECT)
	./$(SELECT) $(JOBID)
runFit: $(FIT)
	./$(FIT) $(JOBID)
runFitERR: $(FITERROR)
	./$(FITERROR) $(JOBID)
runResolution: $(RESOLUTION)
	./$(RESOLUTION) $(JOBID)
runAll: $(FIT)
	./$(MAKEPAIRS)  $(JOBID)
	./$(SELECT) $(JOBID)
	./$(FIT) $(JOBID)


cleanSelect:
	rm -f $(SELECT)

cleanPairs: 
	rm -f $(MAKEPAIRS)
cleanFit:
	rm -f $(FIT)
clean:
	rm -f $(OBJS) core $(FIT) $(SELECT) $(MAKEPAIRS)


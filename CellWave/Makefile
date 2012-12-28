# 
# This Makefile was created by the configure script using
# the files config/MakeDefs.linux and primer/Makefile.in 
# 
# linux with dynamic libraries (requires ecgs compiler)
CELLWAVE = .

SUNDIALS_ARCH = `uname -s`.`uname -m`
#SUNDIALS_LIB = -L$(SUNDIALS)/lib  -lcvode.$(SUNDIALS_ARCH) \
#	      -lshared.$(SUNDIALS_ARCH) -lnvecserial.$(SUNDIALS_ARCH)
#SUNDIALS_INC = -I$(SUNDIALS)/include

CVODE_INC = -I$(SUNDIALS)/include
CVODE_LIB =  -L$(SUNDIALS)/lib  -lcvode.$(SUNDIALS_ARCH) \
         -lshared.$(SUNDIALS_ARCH) -lnvecserial.$(SUNDIALS_ARCH) -lm

#PATCHES= $(Overture)/GridFunction/InterpolatePoints.o
PATCHES= InterpolatePoints.o

#..platform dependent code
include ./Makefile.inc.intel
#include ./Makefile.inc

CWINCLUDES = -I$(CELLWAVE) -I$(BOOST) -IAlpha_cw 

.SUFFIXES:
.SUFFIXES: .C .C .o .f .cc .o

.C.o :;	  $(CC)  $(CWINCLUDES) $(CVODE_INC) $(CCFLAGS) -c $*.C 
.cc.o :;	  $(CC)  $(CWINCLUDES) $(CCFLAGS) -c $*.cc
#-o $*.o
#.cc.o :;  $(CC) -g -I.. $(FTKINCLUDES) -c $*.cc
#.C: $(LIB_DEPENDENCIES) $(PATCHES) #interpolateFluxBoundary.o
	 $(CC) $(CWINCLUDES) $(CCFLAGS) -o $@ $< $(PATCHES) $(CLIBS) $(FLIBS)  $(GLIBS)
.C: 
	@echo "Make Error: you specified a non-existing target."
	@echo "please add a Makefile target for the file you are trying to compile"

## ..CellWave 
CELLWAVERXN_C = GenericReaction.C GenericCalciumIP3Reaction.C \
	         ReactionLiRinzelWagner.C ReactionSlepchenko2Buffer.C Reaction2Buffer.C
CELLWAVESOLVERS_C = GenericSolver.C SolverLiRinzel.C SolverSlepchenko2Buffer.C Solver2Buffer.C
CELLWAVEDRIVER_C = cellWaveMain.C
CELLWAVEOBJ_C = ReactionFactory.C SolverFactory.C getDiffusionDT.C \
	      Info.C DPrintf.C CellWave.C Nucleus.C NucleusGridFunction.C \
	      Probes.C FluxBCBase.C FluxBC.C

CELLWAVE_ODE_C = odeLiRinzel.C odeSlepchenko2Buffer.C

CELLWAVERXN     = ${CELLWAVERXN_C:.C=.o}
CELLWAVESOLVERS = ${CELLWAVESOLVERS_C:.C=.o}
CELLWAVEOBJ     = ${CELLWAVEOBJ_C:.C=.o} $(CELLWAVERXN) $(CELLWAVESOLVERS)
CELLWAVEDRIVER_OBJ = ${CELLWAVEDRIVER_C:.C=.o}

CELLWAVEALL_C    = ${CELLWAVERXN_C} ${CELLWAVESOLVERS_C} ${CELLWAVEOBJ_C} \
                   ${CELLWAVE_ODE_C} ${CELLWAVEDRIVER_C}
CELLWAVE_HEADERS = CellWave.h  GenericReaction.h  getDiffusionDT.h \
	           ReactionFactory.h    SolverLiRinzel.h SolverSlepchenko2Buffer.h \
		   Solver2Buffer.h \
                   DPrintf.h \
                   GenericReactionMacros.h  Info.h \
	           ReactionLiRinzelWagner.h ExplicitSolver.h \
	           GenericSolver.h          ParameterReader.h \
	           ReactionSlepchenko2Buffer.h Reaction2Buffer.h \
	           Probes.h Nucleus.h Info.h FluxBC.h FluxBCBase.h

all = cellwave

cellwave= $(CELLWAVEOBJ) $(CELLWAVEDRIVER_OBJ) $(PATCHES)
cellwave: $(cellwave)
	$(CC) $(CCFLAGS) -o $@ $(cellwave) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

cellWaveMain:
	echo "cannot build cellWaveMain -- try 'make cellwave' instead"

plotCellWave= plotCellWave.o plotCellWaveDriver.o
plotCellWave: $(plotCellWave)
	$(CC) -o plotCellWave $(plotCellWave) \
	$(CLIBS) $(FLIBS) $(GLIBS)

cw_print=$(CELLWAVE_HEADERS) $(CELLWAVEALL_C)
print:
	enscript -2rG -Ecpp $(cw_print) -p - | lp -d cascptr4

cw_source=$(CELLWAVE_HEADERS) $(CELLWAVEALL_C) Makefile *.m
release_name=cw_2003_05_30
release:
	/bin/rm -rf $(release_name)
	mkdir $(release_name)
	mkdir $(release_name)/OUT
	cp $(cw_source) $(release_name)
	cp -r config $(release_name)/config
	cp -r Grids  $(release_name)/Grids
	cp -r PARAMS $(release_name)/PARAMS
	/bin/rm -rf $(release_name)/*/CVS
	tar cvf release$(release_name).tar $(release_name)


ODESOLVEROBJ = DPrintf.o CellWave.o Info.o GenericReaction.o GenericCalciumIP3Reaction.o

odeLiRinzel= odeLiRinzel.o ReactionLiRinzelWagner.o $(ODESOLVEROBJ) 
odeLiRinzel: $(odeLiRinzel)
	$(CC) $(CCFLAGS) -o $@ $(odeLiRinzel) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

odeSlepchenko2Buffer= odeSlepchenko2Buffer.o ReactionSlepchenko2Buffer.o $(ODESOLVEROBJ) 
odeSlepchenko2Buffer: $(odeSlepchenko2Buffer)
	$(CC) $(CCFLAGS) -o $@ $(odeSlepchenko2Buffer) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

radialSolverSlepchenko2Buffer= radialSolverSlepchenko2Buffer.o $(CELLWAVEOBJ)
radialSolverSlepchenko2Buffer: $(radialSolverSlepchenko2Buffer)
	$(CC) $(CCFLAGS) -o $@ $(radialSolverSlepchenko2Buffer) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

ode2Buffer= ode2Buffer.o Reaction2Buffer.o $(ODESOLVEROBJ) 
ode2Buffer: $(ode2Buffer)
	$(CC) $(CCFLAGS) -o $@ $(ode2Buffer) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

radialSolver2Buffer= radialSolver2Buffer.o $(CELLWAVEOBJ)
radialSolver2Buffer: $(radialSolver2Buffer)
	$(CC) $(CCFLAGS) -o $@ $(radialSolver2Buffer) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

matlabParametersFrom2Buffer=matlabParametersFrom2Buffer.o Reaction2Buffer.o $(ODESOLVEROBJ)
matlabParametersFrom2Buffer: $(matlabParametersFrom2Buffer)
	$(CC) $(CCFLAGS) -o $@ $(matlabParametersFrom2Buffer) \
	$(CLIBS) $(FLIBS) $(GLIBS)

#..tests

testFluxBC = testFluxBC.o FluxBC.o FluxBCBase.o getDiffusionDT.o
testFluxBC: $(testFluxBC)
	$(CC) $(CCFLAGS) -o $@ $(testFluxBC) \
	$(CLIBS) $(FLIBS) $(GLIBS)


testInterpGridArray=testInterpGridArray.o InterpGridArray.o
testInterpGridArray: $(testInterpGridArray)
	$(CC) $(CCFLAGS) -o $@ $(testInterpGridArray) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testInterpolatePoints=testInterpolatePoints.o InterpGridArray.o
testInterpolatePoints: $(testInterpolatePoints)
	$(CC) $(CCFLAGS) -o $@ $(testInterpolatePoints) $(PATCHES)\
	$(CLIBS) $(FLIBS) $(GLIBS)

testSubset=testSubset.o
testSubset: $(testSubset)
	$(CC) $(CCFLAGS) -o $@ $(testSubset) $(PATCHES)\
	$(CLIBS) $(FLIBS) $(GLIBS)

testInterpolateMultipleGrids=testInterpolateMultipleGrids.o 
testInterpolateMultipleGrids: $(testInterpolateMultipleGrids)
	$(CC) $(CCFLAGS) -o $@ $(testInterpolateMultipleGrids) $(PATCHES)\
	$(CLIBS) $(FLIBS) $(GLIBS)

testArraySimple= testArraySimple.o
testArraySimple: $(testArraySimple)
	$(CC) $(CCFLAGS) -o $@ $(testArraySimple) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testDiffusion=testDiffusion.o DPrintf.o CellWave.o
testDiffusion: $(testDiffusion)
	$(CC) $(CCFLAGS) -o $@ $(testDiffusion) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testMatchedBoundaryFluxBC=testMatchedBoundaryFluxBC.o DPrintf.o CellWave.o
testMatchedBoundaryFluxBC: $(testMatchedBoundaryFluxBC)
	$(CC) $(CCFLAGS) -o $@ $(testMatchedBoundaryFluxBC) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testDPrintf=testDPrintf.o DPrintf.o CellWave.o
testDPrintf: $(testDPrintf)
	$(CC) $(CCFLAGS) -o $@ $(testDPrintf) \
	$(CLIBS) $(FLIBS) $(GLIBS)

#..CVODE tests
callingCVode= callingCVode.o CVodeOverture.o
callingCVode: $(callingCVode)
	$(CC) $(CCFLAGS) -o $@ $(callingCVode) \
	$(CLIBS) $(FLIBS) $(GLIBS)

cvkx_x= cvkx.c
cvkx.x: $(cvkx_x)
	gcc $(CCFLAGS) $(CVODE_INC) -o cvkx.x $(cvkx_x) \
	$(CVODE_LIB)

cvkx= cvkx.c
cvkx: $(cvkx)
	$(CC) $(CCFLAGS) $(CVODE_INC) -o cvkx $(cvkx) \
	$(CLIBS) $(FLIBS) $(GLIBS) $(CVODE_LIB)

cvkxOverture= cvkxOverture.C
cvkxOverture: $(cvkxOverture)
	$(CC) $(CCFLAGS) $(CVODE_INC) -o cvkxOverture $(cvkxOverture) \
	$(CLIBS) $(FLIBS) $(GLIBS) $(CVODE_LIB)

## ..Flux bc tests
testInterpolateFluxBC     = testInterpolateFluxBC.o getDiffusionDT.o FluxBC.o
testInterpolateFluxBC_dep = $(testInterpolateFluxBC) FluxBC.h getDiffusionDT.h 
testInterpolateFluxBC: $(testInterpolateFluxBC_dep)
	$(CC) $(CCFLAGS) -o testInterpolateFluxBC $(testInterpolateFluxBC) \
	$(CLIBS) $(FLIBS) $(GLIBS) 

testNewFluxBC     = testNewFluxBC.o getDiffusionDT.o FluxBC.o FluxBCBase.o
testNewFluxBC_dep = $(testNewFluxBC) FluxBC.h getDiffusionDT.h 
testNewFluxBC: $(testNewFluxBC_dep)
	$(CC) $(CCFLAGS) -o testNewFluxBC $(testNewFluxBC) \
	$(CLIBS) $(FLIBS) $(GLIBS) 

testNewFluxBC2     = testNewFluxBC2.o getDiffusionDT.o FluxBC.o FluxBCBase.o $(PATCHES)
testNewFluxBC2_dep = $(testNewFluxBC2) FluxBC.h getDiffusionDT.h 
testNewFluxBC2: $(testNewFluxBC2_dep)
	$(CC) $(CCFLAGS) -o testNewFluxBC2 $(testNewFluxBC2) \
	$(CLIBS) $(FLIBS) $(GLIBS) 

testNewFluxBC3     = testNewFluxBC3.o getDiffusionDT.o FluxBC.o FluxBCBase.o $(PATCHES)
testNewFluxBC3_dep = $(testNewFluxBC3) FluxBC.h getDiffusionDT.h 
testNewFluxBC3: $(testNewFluxBC3_dep)
	$(CC) $(CCFLAGS) -o testNewFluxBC3 $(testNewFluxBC3) \
	$(CLIBS) $(FLIBS) $(GLIBS) 

testNewFluxBC4     = testNewFluxBC4.o getDiffusionDT.o FluxBC.o FluxBCBase.o $(PATCHES)
testNewFluxBC4_dep = $(testNewFluxBC4) FluxBC.h getDiffusionDT.h 
testNewFluxBC4: $(testNewFluxBC4_dep)
	$(CC) $(CCFLAGS) -o testNewFluxBC4 $(testNewFluxBC4) \
	$(CLIBS) $(FLIBS) $(GLIBS) 

testNewFluxBC5     = testNewFluxBC5.o getDiffusionDT.o FluxBC.o FluxBCBase.o $(PATCHES)
testNewFluxBC5_dep = $(testNewFluxBC5) FluxBC.h getDiffusionDT.h 
testNewFluxBC5: $(testNewFluxBC5_dep)
	$(CC) $(CCFLAGS) -o testNewFluxBC5 $(testNewFluxBC5) \
	$(CLIBS) $(FLIBS) $(GLIBS) 


## ..OXCWave tests
testReactions:
	echo "Executables have an extension .x"

odeReactionSolver.x: odeReactionSolver.o Reactions.o
	$(CC) -g -o $@ odeReactionSolver.o Reactions.o  -lm

oneDimReactionSolver.x: oneDimReactionSolver.o Reactions.o
	$(CC) -o $@ oneDimReactionSolver.o Reactions.o  -lm



##
## .. alpha_cw
##
proto_20021023= Alpha_cw/proto_20021023.o $(CELLWAVEOBJ)
proto_20021023: $(proto_20021023) 
	$(CC) $(CCFLAGS) -o $@ $(proto_20021023) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

test_20021023 = test_20021023.o  rxnTest.o
test_20021023: $(test_20021023) 
	$(CC) $(CCFLAGS) -o $@ $(test_20021023) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

testGetDiffusionDT = testGetDiffusionDT.o getDiffusionDT.o
testGetDiffusionDT: $(testGetDiffusionDT) 
	$(CC) $(CCFLAGS) -o $@ $(testGetDiffusionDT)  \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

NucleusFiles=Nucleus.o NucleusGridFunction.o

testNucleus= testNucleus.o $(NucleusFiles)
testNucleus: $(testNucleus)
	$(CC) $(CCFLAGS) -o $@ $(testNucleus) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testNucleusGridFunction1= testNucleusGridFunction1.o $(NucleusFiles)
testNucleusGridFunction1: $(testNucleusGridFunction1)
	$(CC) $(CCFLAGS) -o $@ $(testNucleusGridFunction1) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testNucleusGridFunction2= testNucleusGridFunction2.o $(NucleusFiles)
testNucleusGridFunction2: $(testNucleusGridFunction2)
	$(CC) $(CCFLAGS) -o $@ $(testNucleusGridFunction2) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testNucleusGridFunction3= testNucleusGridFunction3.o $(CELLWAVEOBJ)
testNucleusGridFunction3: $(testNucleusGridFunction3)
	$(CC) $(CCFLAGS) -o $@ $(testNucleusGridFunction3) \
	$(CLIBS) $(FLIBS) $(GLIBS)


testProbes= testProbes.o  $(CELLWAVEOBJ)
testProbes: $(testProbes)
	$(CC) $(CCFLAGS) -o $@ $(testProbes) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testProbes2= testProbes2.o $(CELLWAVEOBJ)
testProbes2: $(testProbes2)
	$(CC) $(CCFLAGS) -o $@ $(testProbes2) \
	$(CLIBS) $(FLIBS) $(GLIBS)

testInterpolate=testInterpolate.o $(PATCHES)
testInterpolate: $(testInterpolate)
	$(CC) $(CCFLAGS) -o $@ $(testInterpolate) \
	$(CLIBS) $(FLIBS) $(GLIBS)

##..pre 2002/10/23
proto_20021021= Alpha_cw/proto_20021021.o \
	                Alpha_cw/proto_chem_20021021.o
proto_20021021: $(proto_20021021) 
	$(CC) $(CCFLAGS) -o $@ $(proto_20021021) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

proto_20021018_v4= Alpha_cw/proto_20021018_v4.o \
	                Alpha_cw/proto_chem_20021018_v4.o
proto_20021018_v4: $(proto_20021018_v4) 
	$(CC) $(CCFLAGS) -o $@ $(proto_20021018_v4) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

Alpha_cw/proto_chem_20021018_v4.o: Alpha_cw/proto_chem_20021018_v4.h  Alpha_cw/proto_chem_20021018_v4.C
Alpha_cw/proto_20021018_v4.o: Alpha_cw/proto_chem_20021018_v4.h       Alpha_cw/proto_20021018_v4.C

proto_20021018_v3= Alpha_cw/proto_20021018_v3.o \
	                Alpha_cw/proto_chem_20021018_v3.o
proto_20021018_v3: $(proto_20021018_v3)
	$(CC) $(CCFLAGS) -o $@ $(proto_20021018_v3) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

alpha_cw_2002_10_18_v2= Alpha_cw/alpha_cw_2002_10_18_v2.o \
	             Alpha_cw/alpha_cw_chem_2002_10_18_v2.o
alpha_cw_2002_10_18_v2: $(alpha_cw_2002_10_18_v2)
	$(CC) $(CCFLAGS) -o $@ $(alpha_cw_2002_10_18_v2) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 

alpha_cw_2002_10_18= Alpha_cw/alpha_cw_2002_10_18.o \
	             Alpha_cw/alpha_cw_chem_2002_10_18.o
alpha_cw_2002_10_18: $(alpha_cw_2002_10_18)
	$(CC) $(CCFLAGS) -o $@ $(alpha_cw_2002_10_18) \
	 $(CLIBS) $(FLIBS) $(GLIBS) 



example6-influx = example6-influx.o
example6-influx: $(example6-influx) 
	$(CC) $(CCFLAGS) -o example6-influx $(example6-influx) $(CLIBS) $(FLIBS) $(GLIBS) 


mappedGridExample4 = mappedGridExample4.o ChannelMapping.o
mappedGridExample4: $(mappedGridExample4) 
	$(CC) $(CCFLAGS) -o mappedGridExample4 $(mappedGridExample4) $(CLIBS) $(FLIBS) $(GLIBS) 

mappedGridExample6 = mappedGridExample6.o getDt.o
mappedGridExample6: $(mappedGridExample6) 
	$(CC) $(CCFLAGS) -o mappedGridExample6 $(mappedGridExample6) $(CLIBS) $(FLIBS) $(GLIBS) 

# 2d wave equation
wave = wave.o getDt.o
wave: $(wave) 
	$(CC) $(CCFLAGS) -o wave $(wave) $(CLIBS) $(FLIBS) $(GLIBS) 

# hyperbolic amr example
amrHype= amrHype.o getDt.o 
amrHype: $(amrHype) 
	$(PURIFY) $(CC) $(CCFLAGS) $(PURIFY_OPTIONS) -o amrHype $(amrHype) $(CLIBS) $(FLIBS) $(GLIBS)


clean:  
	rm -f $(PROGRAMS) *.o Alpha_cw/*.o

#depend_date: depend; perl $(Overture)/config/makeDepend -f depend -Y -- $(CC_INCLUDES_FOR_DEPEND) -- $(Source)
#	touch $@
# include the dependencies here
#include depend

#.PRECIOUS: depend_date


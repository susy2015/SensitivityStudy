IDIR       = .
ODIR       = obj
SDIR       = .

CXX        = g++

CXXFLAGS += -I. -std=c++0x -I$(CMSSW_BASE)/src/ -I$(LHAPDF_DATA_PATH)/../../include -I$(CMSSW_RELEASE_BASE)/src 
CXXFLAGS += -I$(shell scram tool info fastjet | grep '^INCLUDE' | sed 's/INCLUDE=//')
CXXFLAGS += -I$(shell scram tool info boost | grep '^INCLUDE' | sed 's/INCLUDE=//')
## Optimization flag
CXXFLAGS += -g -fPIC
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

OPENCV_DIRECTORY=$(CMSSW_BASE)/src/opencv
LIBS       = $(shell root-config --glibs)
LIBS      += -Lobj -lSusyAnaTools
LIBS      += -L$(CMSSW_BASE)/src/TopTagger/TopTagger/test -lTopTagger
LIBS      += -L$(OPENCV_DIRECTORY)/lib/ -lopencv_ml -lopencv_core
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -L$(CMSSW_RELEASE_BASE)/lib/${SCRAM_ARCH}/ -lPhysicsToolsHeppy
LHAPDFLIB  = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF
FILES = $(wildcard $(SDIR)/*.cc)
OBJS  := $(FILES:$(SDIR)/%.cc=$(ODIR)/%.o)
PROGRAMS = tupleTest SSTrimAndSlim_LLHadTau SSTrimAndSlim_ZinvQCD SSTrimAndSlim_TTZ SSTrimAndSlim_Signal

all: mkobj $(OBJS) SusyAnaTools sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

SusyAnaTools : $(OBJS)
	ar rcs lib$@.a $^
	gcc -shared -o $(ODIR)/lib$@.so $^

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

tupleTest: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/PDFUncertainty.o $(ODIR)/tupleReadTest.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

SSTrimAndSlim_LLHadTau: $(ODIR)/SSTrimAndSlim_LLHadTau.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/SSTrimAndSlimCommon.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_ZinvQCD: $(ODIR)/SSTrimAndSlim_ZinvQCD.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/SSTrimAndSlimCommon.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_TTZ: $(ODIR)/SSTrimAndSlim_TTZ.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/SSTrimAndSlimCommon.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_Signal: $(ODIR)/SSTrimAndSlim_Signal.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/SSTrimAndSlimCommon.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(ODIR)/*.so $(PROGRAMS) core 

-include $(ODIR)/*.d

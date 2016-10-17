IDIR       = .
ODIR       = obj
SDIR       = ../../SusyAnaTools/Tools

#CXX       = ~/.vim/bundle/ClangComplete/bin/cc_args.py g++
CXX        = g++

#CXXFLAGS  += -Iinc -std=c++0x -I${CMSSW_BASE}/src/
CXXFLAGS  += -I. -std=c++0x -I$(CMSSW_BASE)/src/ -I$(LHAPDF_DATA_PATH)/../../include

## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD       = g++
LDFLAGS  =

LIBS     = $(shell root-config --glibs)
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2
LHAPDFLIB  = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF

PROGRAMS = tupleTest SSTrimAndSlim_LLHadTau SSTrimAndSlim_ZinvQCD SSTrimAndSlim_TTZ SSTrimAndSlim_Signal

all: mkobj sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

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

$(ODIR)/%.o : $(TTSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TTSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

tupleTest: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/PDFUncertainty.o $(ODIR)/tupleReadTest.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

SSTrimAndSlim_LLHadTau: $(ODIR)/SSTrimAndSlim_LLHadTau.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_ZinvQCD: $(ODIR)/SSTrimAndSlim_ZinvQCD.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_TTZ: $(ODIR)/SSTrimAndSlim_TTZ.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

SSTrimAndSlim_Signal: $(ODIR)/SSTrimAndSlim_Signal.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(ODIR)/*.so $(PROGRAMS) core 

-include $(ODIR)/*.d

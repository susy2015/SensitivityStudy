IDIR       = .
ODIR       = obj
SDIR       = $(CMSSW_BASE)/src/SusyAnaTools/Tools

CXX        = g++
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

PROGRAMS = SS SSDataCardCompare

all: mkobj $(PROGRAMS)

mkobj:
	@mkdir -p obj

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

SS: $(ODIR)/SS.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/searchBins.o $(ODIR)/SSReWeighting.o $(ODIR)/SSDataCard.o $(ODIR)/SignalDataCard.o
	$(LD) $^ $(LIBS)  -o $@

SSDataCardCompare: $(ODIR)/SSDataCardCompare.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS)  -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(ODIR)/*.so $(PROGRAMS) core 

-include $(ODIR)/*.d

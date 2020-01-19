DIR := ${CURDIR}

HDR = ./interface/
SRC = ./src/
PLG = ./plugins/
PRG = ./test/
OBJ = ./obj/
LIB = ./lib/
BIN = ./bin/

HDRSuf = .h
SRCSuf = .cc
PRGSuf = .cpp
OBJSuf = .o
LIBSuf = .so
BINSuf = .exe

HDRS    =  $(wildcard $(HDR)*$(HDRSuf))
SRCS    =  $(wildcard $(SRC)*$(SRCSuf))
PLGS    =  $(wildcard $(PLG)*$(HDRSuf))
_OBJS   =  $(patsubst %$(SRCSuf), %$(OBJSuf), $(SRCS))
OBJS    =  $(patsubst $(SRC)%, $(OBJ)%, $(_OBJS))
PRGS    =  $(wildcard $(PRG)*$(PRGSuf))
_BINS   =  $(wildcard $(PRG)*$(PRGSuf))
__BINS  =  $(_BINS:$(PRGSuf)=$(BINSuf))
___BINS =  $(notdir $(__BINS))
BINS    =  $(addprefix $(BIN),${___BINS})

LINKDEF   =  $(wildcard ${HDR}*LinkDef.h)
DICTHDRS  =  $(patsubst $(LINKDEF),,$(HDRS)) $(LINKDEF)


ARCH  =  $(shell root-config --arch)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs) -lGenVector -lFoam -lMinuit -lTMVA -lMLP -lXMLIO  -lTreePlayer



CXX  =  g++
CXXFLAGS  = -Wall -O2 -fPIC -I$(DIR) $(ROOTCFLAGS) 

CPP  =  g++
CPPFLAGS  = -Wall -I$(DIR) $(ROOTCFLAGS)

LD       =  g++
LDFLAGS  =  -rdynamic -shared -O2
SONAME	 =  libCfgManager.so
SOFLAGS  =  -Wl,-soname,

GLIBS   =  -lm -ldl -rdynamic $(ROOTGLIBS)



#################################################
#if mac 64
ifeq ($(ARCH),macosx64)
LIBSuf  =  .dylib

CPPFLAGS  =  -Wall -W -Woverloaded-virtual -O2 -pipe -I$(DIR) $(ROOTCFLAGS)

CXXFLAGS  =  -Wall -W -Woverloaded-virtual -O2 -pipe -I$(DIR) $(ROOTCFLAGS)

LDFLAGS  =  -dynamiclib -shared -single_module -undefined dynamic_lookup
SONAME	 =  libCfgManager.dylib
SOFLAGS  =  -Wl,-install_name,
endif
#################################################



.PHONY: all clean test


all: $(LIB)$(SONAME)

exe: $(BINS)

test:
	@echo "HDRS = $(HDRS)"
	@echo "DICTHDRS = $(DICTHDRS)"
	@echo "SRCS = $(SRCS)"
	@echo "PLGS = $(PLGS)"
	@echo "PRGS = $(PRGS)"
	@echo "OBJS = $(OBJS)"
	@echo "BINS = $(BINS)"

$(BIN)%$(BINSuf): $(PRG)%$(PRGSuf) $(HDRS) $(LIB)$(SONAME)
	$(CPP) $(CPPFLAGS) $(GLIBS) -L$(LIB) -lCfgManager -o $@ $<

$(OBJ)%$(OBJSuf): $(SRC)%$(SRCSuf)
	$(CXX) -c $(CXXFLAGS) -o $@ $< 

$(LIB)mydict.cc: $(DICTHDRS)
	@echo "Generating dictionary..."
	rootcling -f $(LIB)mydict.cc -c -p ${CXXFLAGS} $(DICTHDRS)

$(LIB)mydict.o: $(LIB)mydict.cc 
	$(CXX) -c $(CXXFLAGS) -o $@ $<

$(LIB)$(SONAME): $(OBJS) $(LIB)mydict.o
	@echo "Linking $(SONAME):"
	$(LD) $(LDFLAGS) $(OBJS) $(LIB)mydict.o -o $(LIB)$(SONAME) $(SOFLAGS)$(SONAME)

clean:
	@echo "cleaning..."
	rm -f $(OBJ)*$(OBJSuf) $(LIB)*$(LIBSuf) $(LIB)mydict* $(BIN)*$(BINSuf)

# From command line: export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
MU_PARSERX_PATH = ./external/muparserx/build

# Compiler
CXX = g++

# Optimize flags
OPTFLAGS = -O3
GETFEM_PATH = $(HOME)/getfem-5.4
INCLUDE = -I$(GETFEM_PATH)/include -I$(GETFEM_PATH)/src -I$(GETFEM_PATH)/src/gmm -I./include -I/usr/include -I/usr/include/x86_64-linux-gnu -I./external/muparserx/parser
CXXFLAGS = $(INCLUDE) $(OPTFLAGS) 

# Executable source
EXESRCS = main.cpp

# Executable object file
EXEOBJS = $(EXESRCS:.cpp = .o)

# Executable name
EXEC = main

# Sources folder
FOLDER = src/


# Laptop
LIB_PATH = /usr/lib /usr/lib/x86_64-linux-gnu $(MU_PARSERX_PATH)


# Laptop
LDLIBS = /usr/local/lib/libgetfem.a
####### ADDED
LDLIBS += -Wl,-rpath,/usr/lib/x86_64-linux-gnu
LDLIBS += -L$(MU_PARSERX_PATH) -lmuparserx

######## END ADDED

# Laptop
LDLIBS += $(GETFEM_LIB) -rdynamic /usr/lib/x86_64-linux-gnu/libqhull.so.8.0 /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/libblas.so.3


DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS

# Sources
SRCS = $(wildcard $(FOLDER)*.cpp)

# Objects
OBJS = $(SRCS:.cpp=.o)

# Headers
HEADERS = $(SRCS:.cpp=.hpp)

# Name file of dependences
DEPEND = make.dep

.PHONY: all clean

all : $(DEPEND) $(OBJS) $(EXEOBJS)
	$(CXX) $(OPTFLAGS) -o $(EXEC) $(EXEOBJS) $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE) 

mecc :  $(DEPEND) $(OBJS) main_mecc.o
	$(CXX) $(OPTFLAGS) -o $(EXEC) main_mecc.o $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

$(DEPEND) : $(SRCS) $(EXESRCS)
	$(CXX) -MM $(SRCS) $(EXESRCS) -MF $(DEPEND)  $(INCLUDE) 

-include $(DEPEND)

clean:
	@$(RM) $(OBJS) $(EXECOBJS)

distclean: clean
	$(RM) $(EXEC)
	$(RM) -f ./doc
	$(RM) *.out *.bak *.dep *.log *~ 

doc:
	doxygen $(DOXYFILE)


# Directories
SRC_DIR = ./src
INC_DIR = ./include
BUILD_DIR = ./build
MU_PARSERX_DIR = ./external/muparserx/build

# Set compiler variables
CXX = mpic++
CPPFLAGS = -I$(INC_DIR) \
  -I/usr/local/include \
  -I/usr/local/include/getfem \
  -I/usr/local/include/gmm \
  -I./external/muparserx/parser \
  -I/usr/include \
  -I/usr/include/x86_64-linux-gnu
DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS -DGMM_USES_MPI=0
CXXFLAGS = -std=c++20 -O3 $(CPPFLAGS) $(DEFTAGS) -MMD -MP

# Linker flags: library search paths + runtime linker path
LDFLAGS = \
  -L/usr/lib \
  -L/usr/lib/x86_64-linux-gnu \
  -L$(MU_PARSERX_DIR) \
  -Wl,-rpath,/usr/lib/x86_64-linux-gnu

# Libraries to link against
LDLIBS = \
  -lgetfem \
  -lmuparserx \
  -ldmumps_seq -lmumps_common_seq -lsmumps_seq -lcmumps_seq \
  -llapack -lblas \
  -lqhull \
  -rdynamic


# Get all source files in the src directory
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

# Executable source
EXE_SRCS = main.cpp
EXE_OBJS = $(EXE_SRCS:.cpp=.o)
EXEC = $(EXE_SRCS:.cpp=)

# Dependencies
DEPS = $(OBJS:.o=.d) $(EXE_OBJS:.o=.d)

.PHONY: all clean
.DEFAULT_GOAL = all

all: $(EXEC)

$(EXEC): $(EXE_OBJS) $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	@$(RM) -f $(OBJS) $(EXE_OBJS) $(DEPS)

distclean: clean
	$(RM) $(EXEC)
	$(RM) -f ./doc
	$(RM) *.out *.bak *.log *~ *.vtk

doc:
	doxygen $(DOXYFILE)


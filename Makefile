# make sure there aren't trailing spaces for any directory names! 

# user specific variables 
# root directory
HOME = .
# superlu stuff. where is it installed and what is the library called 
SUPERLU = /opt/superlu
SUPERLULIB = superlu_5.1
# compiler 
CXX = g++

# directories with source files in them 
UTILS = $(HOME)/utils
SRC = $(HOME)/src
FE = $(SRC)/fem
OPERATORS = $(SRC)/operators
SOLVERS = $(SRC)/solvers
TIMER = $(UTILS)/timer
VISIT = $(UTILS)/VisitWriter
EIGEN = $(UTILS)/eigen
MATERIALS = $(SRC)/materials
FIELD = $(SRC)/field

# where make searches for source files 
VPATH = $(FE) $(OPERATORS) $(SOLVERS) $(TIMER) \
	$(VISIT) $(EIGEN)/src $(MATERIALS) $(FIELD)

# where to look for header files 
CFLAGS = -std=c++14 -I$(FE) -I$(OPERATORS) \
	-I$(SOLVERS) -I$(TIMER) -I$(VISIT) \
	-I$(EIGEN) -I$(MATERIALS) -I$(FIELD)
ifdef SUPERLU
CFLAGS += -DSUPERLU
endif

# optimzation or debug 
CFLAGS += -O3
# CFLAGS += -g -Wall

# store object files and dependency files 
OBJ = $(HOME)/obj
DEP = $(HOME)/dep

# libraries we are dependent on: SuperLU (which depends on BLAS) 
ifdef SUPERLU
LIBS = -I$(SUPERLU)/SRC -L$(SUPERLU)/lib -l$(SUPERLULIB)
LIBS += -lblas
endif 

# get all file names for all .cpp files 
SRCFILES = $(notdir $(wildcard $(FE)/*.cpp $(OPERATORS)/*.cpp \
	$(SOLVERS)/*.cpp $(TIMER)/*.cpp $(VISIT)/*.cpp\
	$(MATERIALS)/*.cpp $(FIELD)/*.cpp))

# convert to object files and dependency files 
OBJS = $(patsubst %.cpp,$(OBJ)/%.o,$(SRCFILES))
DEPS = $(patsubst $(OBJ)/%.o,$(DEP)/%.d, $(OBJS))

$(OBJ)/%.o : %.cpp $(HOME)/Makefile 
	mkdir -p $(OBJ); $(CXX) -c $(CFLAGS) $(LIBS) $< -o $@ 
	mkdir -p $(DEP)
	$(CXX) -MM $(CFLAGS) $(LIBS) $< | sed -e '1s@^@$(OBJ)\/@' > $*.d; mv $*.d $(DEP)

-include $(DEPS)

cleanall : 
	rm -rf $(OBJ); rm -rf $(DEP); rm -f *.exe; rm -f *.vtk; rm -f time.table

clean : 
	rm -f *.exe *.vtk time.table residual err matrix 

listsrc : 
	@echo $(SRCFILES) 
listobj : 
	@echo $(OBJS)
listdep :
	@echo $(DEPS) 
listflags : 
	@echo $(CFLAGS)

.PHONY : docs
docs : 
	cd $(HOME)/docs; doxygen Doxyfile

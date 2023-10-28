CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp

DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_INCLUDE = ./include
DIR_DATA = ./data

SRC = $(wildcard $(DIR_SRC)/*.cpp)
OBJS = $(patsubst %.cpp,$(DIR_OBJ)/%.o,$(notdir $(SRC)))
INC = $(wildcard $(DIR_INCLUDE)/*.hpp)

all: solver

solver: $(OBJS)
	$(CXX) $(CXXFLAGS) -o solver $(OBJS)

$(OBJS): $(SRC) $(INC)
	@echo $^
	$(CXX) $(CXXFLAGS) -c $(DIR_SRC)/$(notdir $*).cpp -I $(DIR_INCLUDE) -o $@

clean:
	rm -f *.o $(OBJS)
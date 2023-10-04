DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_INCLUDE = ./include

SRC = $(wildcard $(DIR_SRC)/*.cpp)
OBJS = $(patsubst %.cpp,$(DIR_OBJ)/%.o,$(notdir $(SRC)))
INC = $(wildcard $(DIR_INCLUDE)/*.hpp)

all: model
	@echo SRC $(SRC)
	@echo OBJS $(OBJS)
	@echo notdir_SRC $(notdir $(SRC))
	@echo INC $(INC)

model: $(OBJS)
	g++ -o model $(OBJS)

$(OBJS): $(SRC) $(INC)
	@echo $^
	g++ -c $(DIR_SRC)/$(notdir $*).cpp -I $(DIR_INCLUDE) -o $@

clean:
	rm -f *.o $(OBJS)
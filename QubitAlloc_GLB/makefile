CC = g++

CPP_VER = -std=c++17
CFLAG = -Wall -c
OFLAG = -O3

SRC_DIR = src
INC_DIR = include

EXEC = qabb

# OS detection
ifeq ($(OS), Windows_NT) 
    RM = del /Q
    EXEC_EXT = .exe
else
    RM = rm -f
    EXEC_EXT =
endif

EXEC_FILES = $(addsuffix $(EXEC_EXT),$(EXEC))

all: $(EXEC)

## executable
qabb: main.o node.o preprocessing.o hungarian.o objective.o utils.o
	$(CC) $(CPP_VER) $(OFLAG) $^ -o $@

## objects
main.o: $(SRC_DIR)/main.cpp $(INC_DIR)/*.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<

node.o: $(SRC_DIR)/node.cpp $(INC_DIR)/*.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<

preprocessing.o: $(SRC_DIR)/preprocessing.cpp $(INC_DIR)/*.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<

objective.o: $(SRC_DIR)/objective.cpp $(INC_DIR)/*.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<

hungarian.o: $(SRC_DIR)/hungarian.cpp $(INC_DIR)/*.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<

utils.o: $(SRC_DIR)/utils.cpp $(INC_DIR)/utils.hpp
	$(CC) $(CPP_VER) $(OFLAG) $(CFLAG) $<


## clean
.PHONY: clean fclean

clean:
	$(RM) *.o

fclean: clean
	$(RM) $(EXEC_FILES)

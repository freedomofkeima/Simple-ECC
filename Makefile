# Global configuration
CC = gcc
EXT = c
LIB = -lgmp
MODULES = $(wildcard src/*.$(EXT))
FLAGS += -O2 -g

# Compiling main
OBJS := ${MODULES:src/%.$(EXT)=bin/%.o}

all: $(OBJS)
	$(CC) -o bin/main $(OBJS) $(LIB)

# Run
run:
	bin/main

# Cleaning compilation
clean:
	rm -rf bin/*.o bin/*~

bin/%.o : src/%.$(EXT)
	$(CC) -c $< -o $@ $(FLAGS) $(INCLUDE)

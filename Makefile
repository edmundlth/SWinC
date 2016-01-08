CC = gcc
CFLAGS = -Wall
OBJECTS = swnn.o thermodynamics_routines.o scoring_routines.o
SRC = swnn.c thermodynamics_routines.c scoring_routines.c 
DEP = swnn.h
EXE = swnn
EXE_DEBUG = swnn_debug

$(EXE): $(OBJECTS)
		$(CC) $(CFLAGS) -O2 -o $(EXE) $(OBJECTS)

$(EXE_DEBUG): $(DEP)
		$(CC) $(CFLAGS) -g -o $(EXE_DEBUG) $(SRC)

$(OBJECTS):$(DEP)

clean:
		/bin/rm -fr $(EXE) $(EXE_DEBUG) $(OBJECTS) *.dSYM

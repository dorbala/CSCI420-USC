# Makefile 
# CSCI 480
# Assignment 3

# we assume the pic directory locates one level above,
# change PIC_PATH if this is not the case
PIC_PATH = $(abspath $(CURDIR)/../pic)

INCLUDE = -I$(PIC_PATH)
LIBRARIES = -L$(PIC_PATH) -lGL -lGLU -lglut -lpicio -ljpeg -lm

COMPILER = g++
COMPILERFLAGS = -O3 $(INCLUDE)

PROGRAM = ShyamaDorbala_assign3
SOURCE = ShyamaDorbala_assign3.cpp
OBJECT = ShyamaDorbala_assign3.o

.cpp.o: 
	$(COMPILER) -c $(COMPILERFLAGS) $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECT)
	$(COMPILER) $(COMPILERFLAGS) -o $(PROGRAM) $(OBJECT) $(LIBRARIES)

clean:
	-rm -rf core *.o *~ "#"*"#" $(PROGRAM)

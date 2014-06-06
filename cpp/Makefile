CC=mpic++ -std=c++11 -lm -g
CFLAGS=-c -Wall
LFLAGS=
SOURCES=RunnerModule.cpp mpi.cpp simul_constants.cpp
OBJECTS=$(SOURCES:.cpp=.o)
CLEAN= $(OBJECTS)
EXECUTABLE=simul

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(CLEAN)

rebuild: clean all

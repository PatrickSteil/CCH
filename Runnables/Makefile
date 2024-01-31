CC=g++ -fopenmp
FLAGS=-std=c++20 -pipe
OPTIMIZATION=-march=native -O3
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wfatal-errors -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer
RELEASE=-ffast-math -ftree-vectorize -Wfatal-errors -DNDEBUG -flto -fomit-frame-pointer -fno-stack-protector -mtune=native -Ofast
EXTRAFLAGS=-lpthread

all: Release

Release:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o CCHTest CCHTest.cpp $(EXTRAFLAGS)

Debug:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o CCHTest CCHTest.cpp $(EXTRAFLAGS)

ReadGraphRelease:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o ReadGraph ReadGraph.cpp $(EXTRAFLAGS)

ReadGraphDebug:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o ReadGraph ReadGraph.cpp $(EXTRAFLAGS)

ExperimentsRelease:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o Experiments Experiments.cpp $(EXTRAFLAGS)

ExperimentsDebug:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o Experiments Experiments.cpp $(EXTRAFLAGS)

EXE = postsort
OBJ = main.o BinaryTreeSort.o

CXX = g++

INCL = state.h BinaryTreeSort.h

$(EXE): $(OBJ)
	$(CXX) $(OBJ) -o $(EXE)

main.o: main.c++ $(INCL)
	$(CXX) -c main.c++
BinaryTreeSort.o: BinaryTreeSort.cc $(INCL)
	$(CXX) -c BinaryTreeSort.cc

clean:
	rm -rf *.o

install:
	cp -r $(EXE) ../../bin/

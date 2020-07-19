WARN_FLAGS = \
	-std=c++11 -O3 -g  \
	 -I/Users/rodtheo/Bioinfo/PROJETOS/STUDIES/LZ/sdsl-lite/build/include \
	  -I/Users/rodtheo/Bioinfo/PROJETOS/STUDIES/LZ/BWTIL/extern/bitvector/include \
		 -L/Users/rodtheo/Bioinfo/PROJETOS/STUDIES/LZ/sdsl-lite/build/lib -lsdsl -ldivsufsort -ldivsufsort64

BIN = sFM-index-sdsl.x
OBJECTS = sFM-index-sdsl.o

#BIN2 = test2
#OBJECTS2 = test2.o

#all : $(BIN) $(BIN2)
all: $(BIN)

$(BIN): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(WARN_FLAGS)

#$(BIN2): $(OBJECTS2)
	#$(CXX) -o $@ $(OBJECTS2) $(WARN_FLAGS)

sFM-index-sdsl.o: sFM-index-vSDLS.cpp optparse.h
	$(CXX) $(WARN_FLAGS) -c $< -o $@

# Check includes in multiple files.
check: test.o
	echo '#include "optparse.h"' > tmp_second.cc
	$(CXX) -o tmp_test test.o tmp_second.cc
	rm tmp_second.cc
	./tmp_test
	rm tmp_test

clean:
	rm -f *.o $(BIN) $(BIN2)

.PHONY: check clean

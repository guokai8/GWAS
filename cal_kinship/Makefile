CXX = g++
CPPFLAGS = -I/arma/include 
CXXFLAGS = -O2 -Wall -DARMA_DONT_USE_WRAPPER 
LDFLAGS = -lblas -llapack
SRC_DIR=src
kin_cal : $(SRC_DIR)/main.cpp 
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(SRC_DIR)/main.cpp
.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $c
clean:
	rm kin_cal



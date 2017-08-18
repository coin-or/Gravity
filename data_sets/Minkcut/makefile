SYSTEM     = x86-64_darwin
CCC = g++ -std=c++11  -Debug

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCFLAGS = -Wall 

band: band.o 
	$(CCC)  $(CCFLAGS) band.o   -o band 
band.o: band.cpp
	$(CCC)  $(CCFLAGS) -c  band.cpp  -o band.o

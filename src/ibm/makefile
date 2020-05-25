CXX=g++
CXXFLAGS=-Wall -O3 -std=c++11 -ggdb
EXE= xcue_integration.exe
$(EXE) : cue_integration_social_learning.cpp auxiliary.hpp individual.hpp
	$(CXX) $(CXXFLAGS) -o $(EXE) cue_integration_social_learning.cpp individual.cpp


clean :
	rm -rf $(EXE)

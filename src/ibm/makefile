CXX=g++
CXXFLAGS=-Wall -O3 -std=c++11 -ggdb

xcue_integration : cue_integration_social_learning.cpp auxiliary.hpp individual.hpp
	$(CXX) $(CXXFLAGS) -o xcue_integration cue_integration_social_learning.cpp individual.cpp


clean :
	rm -rf xcue_integration

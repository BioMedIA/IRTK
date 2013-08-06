# Makefile

CPP = g++
CPPFLAGS = -w
#CPPFLAGS = -g -Wall
#LDFLAGS = -lboost_thread -lboost_system -lboost_filesystem

maxflow: *.cpp
	$(CPP) $(CPPFLAGS) graph.cpp maxflow.cpp example.cpp -o example $(LDFLAGS)

clean:
	rm -f example

python:
	python setup.py build_ext --inplace

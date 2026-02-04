#FLAGS = -std=c++14 -I`root-config --incdir` -I/home/bpc135/.local/include/ -O3 -g
FLAGS = -std=c++14 -I`root-config --incdir` -I $(HOME)/.local/include/ -O3 -g
ROOTFLAGS = -L`root-config --libdir` `root-config --glibs --new`
#LDFLAGS = $(ROOTFLAGS) -L/home/bpc135/.local/lib/ -lpixie -lClarionTrinity
LDFLAGS = $(ROOTFLAGS) -L $(HOME)/.local/lib/ -lpixie -lClarionTrinity
CXX = g++

all : BasicSort FESort 

BasicSort : BasicSort.cc 
	$(CXX) $(FLAGS) BasicSort.cc -o ../BasicSort $(LDFLAGS)

FESort : FESort.cc 
	$(CXX) $(FLAGS) FESort.cc -o ../FESort $(LDFLAGS)

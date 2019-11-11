objects =getDistributionEncode/DistributionEncode.o 
ARMA_INCLUDE_FLAG = -I   ArmaInclude/include   -I    MlpackInclude  
LIB_FLAGS = -lblas -llapack 
CXXFLAGS = $(ARMA_INCLUDE_FLAG)
CC = g++  -std=c++11  -fopenmp     -g  -I getDistributionEncode/    #-Wall
all:    prepareall
	$(CC)$(CXXFLAGS)   -o   StructureEncode    main.cpp  $(objects) $(LIB_FLAGS)   -lmlpack    
prepareall:    subsystem
subsystem:
	$(MAKE) -C getDistributionEncode/
clean :  cleansub
	rm   StructureEncode
cleansub :
	rm  $(objects)


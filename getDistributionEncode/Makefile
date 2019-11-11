objects = DistributionEncode.o
CC = g++  -std=c++11  -g   #-Wall
ARMA_INCLUDE_FLAG = -I ../ArmaInclude/include  -I ../MlpackInclude
LIB_FLAGS = -lblas -llapack  -lmlpack 
CXXFLAGS = $(ARMA_INCLUDE_FLAG)
install : $(objects)
DistributionEncode.o : DistributionEncode.cpp   DistributionEncode.h   -lmlpack 
	$(CC)$(CXXFLAGS)  -c  $<   $(LIB_FLAGS)
clean :
	rm  $(objects) 

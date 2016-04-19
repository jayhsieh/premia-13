CXX = g++ 
CXXFLAGS+= -Wall
LIBS= -lm 

lmm_bermudaprice_andersen_obj = lmm_bermudaprice_andersen.o\
 mathsb.o\
 lmm_bermuda_andersen_example.o


%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $<


lmm_bermuda_andersen_example : $(lmm_bermudaprice_andersen_obj)
	$(CXX) -o lmm_bermuda_andersen_example $(lmm_bermudaprice_andersen_obj) $(LIBS)

lmm_bermudaprice_andersen.o : lmm_bermudaprice_andersen.cpp lmm_bermudaprice_andersen.h mathsb.h 
lmm_bermuda_andersen_example.o : lmm_bermuda_andersen_example.cpp lmm_bermudaprice_andersen.h mathsb.h 
mathsb.o: mathsb.cpp mathsb.h 


clean :
	rm -f mathsb.o lmm_bermudaprice_andersen.o lmm_bermuda_andersen_example.o

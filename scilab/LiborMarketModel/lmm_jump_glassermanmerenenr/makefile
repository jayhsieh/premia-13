CXX = g++
CXXFLAGS += -Wall
LIBS= -lm 



%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $<


do: lmm_jump_caplet_mc_example lmm_jump_caplet_glassmer_example

lmm_jump_caplet_mc_example: lmm_jump_capletprice_glassmer.o mathsb.o lmm_jump_caplet_mc_example.o
	$(CXX) $(CXXFLAGS) lmm_jump_capletprice_glassmer.o mathsb.o lmm_jump_caplet_mc_example.o -o lmm_jump_caplet_mc_example $(LIBS)

lmm_jump_caplet_glassmer_example: lmm_jump_capletprice_glassmer.o mathsb.o lmm_jump_caplet_glassmer_example.o
	$(CXX) $(CXXFLAGS) lmm_jump_capletprice_glassmer.o mathsb.o lmm_jump_caplet_glassmer_example.o -o lmm_jump_caplet_glassmer_example $(LIBS)



lmm_jump_capletprice_glassmer.o : lmm_jump_capletprice_glassmer.cpp lmm_jump_capletprice_glassmer.h mathsb.h 
lmm_jump_caplet_mc_example.o : lmm_jump_caplet_mc_example.cpp lmm_jump_capletprice_glassmer.h mathsb.h 
lmm_jump_caplet_glassmer_example.o : lmm_jump_caplet_glassmer_example.cpp lmm_jump_capletprice_glassmer.h mathsb.h 
mathsb.o     : mathsb.cpp mathsb.h 



clean :
	rm -f mathsb.o lmm_jump_capletprice_glassmer.o lmm_jump_caplet_mc_example.o lmm_jump_caplet_glassmer_example.o

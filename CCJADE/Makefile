CC = g++
#CFLAGS = -g -m32 -mfpmath=sse -Ofast -flto -march=native -funroll-loops -Wfatal-errors -Wno-deprecated -D__USE_MINGW_ANSI_STDIO
CFLAGS = -g -O3 -Wfatal-errors -Wno-deprecated -D__USE_MINGW_ANSI_STDIO
INCLUDEDIRS = -I./tclap/include -I./qpa -I./rbfn -I./SOCO_SI -I./eigen -I./km -I./libgp/include -I./dlib-18.17
CCDEOBJS = main.o CCDE.o Decomposer.o JADE.o sobol.o QuadraticRegression.o KmeansPP.o RBFNetwork.o

LIBGPOBS = cg.o cov.o cov_factory.o cov_linear_ard.o cov_linear_one.o cov_matern3_iso.o cov_matern5_iso.o cov_noise.o cov_periodic.o cov_periodic_matern3_iso.o cov_prod.o cov_rq_iso.o cov_se_ard.o cov_se_iso.o cov_sum.o gp.o gp_utils.o input_dim_filter.o rprop.o sampleset.o

LIBGPSRC = libgp/src/cg.cc libgp/src/cov.cc libgp/src/cov_factory.cc libgp/src/cov_linear_ard.cc libgp/src/cov_linear_one.cc libgp/src/cov_matern3_iso.cc libgp/src/cov_matern5_iso.cc libgp/src/cov_noise.cc libgp/src/cov_periodic.cc libgp/src/cov_periodic_matern3_iso.cc libgp/src/cov_prod.cc libgp/src/cov_rq_iso.cc libgp/src/cov_se_ard.cc libgp/src/cov_se_iso.cc libgp/src/cov_sum.cc libgp/src/gp.cc libgp/src/gp_utils.cc libgp/src/input_dim_filter.cc libgp/src/rprop.cc libgp/src/sampleset.cc

FUNCTOBJS = funsoft.o

saccjade.exe : $(CCDEOBJS) $(KMOBJS) $(LIBGPOBS) $(FUNCTOBJS)
	$(CC) $(CFLAGS) $(CCDEOBJS) $(KMOBJS) $(LIBGPOBS) $(FUNCTOBJS) -o saccjade.exe

$(CCDEOBJS) : main.cpp CCDE.cpp Decomposer.cpp JADE.cpp sobol.cpp ./dlib/all/source.cpp ./qpa/QuadraticRegression.cpp ./rbfn/KmeansPP.cpp ./rbfn/RBFNetwork.cpp
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c main.cpp CCDE.cpp Decomposer.cpp JADE.cpp sobol.cpp ./dlib/all/source.cpp ./qpa/QuadraticRegression.cpp ./rbfn/KmeansPP.cpp ./rbfn/RBFNetwork.cpp
	
$(FUNCTOBJS) : SOCO_SI/funsoft.cpp
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c SOCO_SI/funsoft.cpp
	
$(LIBGPOBS) : $(LIBGPSRC)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $(LIBGPSRC)
	
.PHONY: clean

clean:
	del *.o *~ core 
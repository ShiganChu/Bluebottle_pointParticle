PREC = DOUBLE

C = mpicc.openmpi
COPT = -std=c99 -pedantic -Wall -Wextra -fopenmp -D$(PREC)
LDINCS = -I /usr/lib/openmpi/include -I /usr/local/cgns/include
LDLIBS = -lm -L /usr/local/hdf5/lib -L /usr/local/cgns/lib -lcgns -lhdf5 

CUDA = nvcc
  # compiled to sm_30 because cuda-memcheck can't handle sm_35
#CUDAOPT = -arch=sm_30 -Xcompiler -fopenmp -m64 -D$(PREC)
CUDAOPT = -arch=sm_30 -Xcompiler -m64 -D$(PREC)
SDK_DIR = /home/asiera/NVIDIA_CUDA-5.0_Samples
CUDA_DIR = /usr/local/cuda/lib64

CUDAINCS = -I $(SDK_DIR)/common/inc -I $(CUDA_DIR)/include -I /home/asiera/NVIDIA_CUDA-5.0_Samples/0_Simple/simplePrintf/
CUDALIBS = -L $(SDK_DIR)/lib -L $(CUDA_DIR) 	\
	-lcudart

SRCC =	bluebottle.c	\
	domain.c	\
	point.c	\
	precursor.c	\
	recorder.c	\
	seeder.c	\
	vtk.c	

SRCCUDA = cuda_bluebottle.cu	\
	cuda_bicgstab.cu	\
	cuda_point.cu	\
	cuda_testing.cu		\
	entrySearch.cu		\
	bluebottle_kernel.cu	\
	bicgstab_kernel.cu	\
	entrySearch_kernel.cu	\
	point_kernel.cu	\
	shigan.cu	\
	shigan_kernel.cu

EXTRA = Makefile		\
	bluebottle.h		\
	cuda_bluebottle.h	\
	cuda_bicgstab.h		\
	cuda_point.h		\
	cuda_testing.h		\
	domain.h		\
	entrySearch.h		\
	point.h		\
	precursor.h		\
	recorder.h		\
	vtk.h		\
	shigan.h

# compile normally
all: COPT += -O3
all: CUDAOPT += -O3
all: bluebottle

# compile for batch job submission
batch: COPT += -O3 -DBATCHRUN
batch: CUDAOPT += -O3
batch: bluebottle

# compile with stair-stepped interior boundaries
steps: COPT += -DSTEPS -O2
steps: CUDAOPT += -DSTEPS -O2
steps: bluebottle

# compile with debug output
#debug: COPT += -DDEBUG -g -ggdb3 -gstabs3 -gdwarf-2
#debug: CUDAOPT += -DDEBUG -g -G   -gencode arch=compute_20,code=sm_20 -v -DTHRUST_DEBUG  

debug: COPT += -DDEBUG -g 
debug: CUDAOPT += -DDEBUG -g -G 
debug: bluebottle

# compile with testing code
test: COPT += -DDEBUG -DTEST -g
test: CUDAOPT += -DDEBUG -DTEST -g -G
test: bluebottle

# write robodoc documentation
doc:
	cd .. && robodoc --html --multidoc --doc doc/robodoc && robodoc --latex --singledoc --sections --doc doc/LaTeX/Bluebottle_0.1_robodoc && cd doc/LaTeX && pdflatex Bluebottle_0.1_robodoc.tex && pdflatex Bluebottle_0.1_robodoc.tex && pdflatex Bluebottle_0.1_robodoc.tex && echo '\nmake doc: Complete.'

OBJS = $(addsuffix .o, $(basename $(SRCC)))
OBJSCUDA = $(addsuffix .o, $(basename $(SRCCUDA)))

%.o:%.cu
	$(CUDA) $(CUDAOPT) -dc $< $(CUDAINCS) $(LDINCS)

%.o:%.c
	$(C) $(COPT) -c $< $(LDINCS)

bblib.o: $(OBJSCUDA)
	$(CUDA) $(CUDAOPT) -dlink $+ -o $@ $(CUDALIBS)

bluebottle: $(OBJSCUDA) bblib.o $(OBJS)
	$(C) $(COPT) -o $@ $+ $(LDLIBS) $(CUDALIBS) -lstdc++

clean:
	rm -f *.o bluebottle seeder

#--------------------------------------------------
# Make file for UMOA calculation
#--------------------------------------------------
COMPILER=GNU
TARGET=test
#--------------------------------------------------
# for compile
#--------------------------------------------------
#LIBRARY_PATH=-L../lib/
INCLUDE_PATH=-I/Users/miyagi/LinAlgLib
FDEP=makedepf90
FC=gfortran
LIBS=-lLinAlg

OMP=-fopenmp
FDFLAGS=-fbounds-check -Wall -fbacktrace -O -Wuninitialized -Ddebug


#--------------------------------------------------
# Source Files
#--------------------------------------------------
SRCDIR = .
DEPDIR = .
SRCF90 += $(wildcard $(SRCDIR)/*.f90)
SRCF95 += $(wildcard $(SRCDIR)/*.F90)
DEPC = $(SRCC:$(SRCDIR/%.c=$(DEPDIR)/%.d))
SRCS = $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95)

OBJDIR = .
OBJF90 += $(SRCF90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJF95 += $(SRCF95:$(SRCDIR)/%.F90=$(OBJDIR)/%.o)
OBJS = $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95)

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: $(TARGET)
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $(TARGET).exe $^
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $< $(INCLUDE_PATH) $(LIBS)

dep:
	$(FDEP) $(SRCS) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

clean:
	rm -f $(TARGET).exe
	rm -f *.mod *.o

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)


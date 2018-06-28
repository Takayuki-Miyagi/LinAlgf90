#--------------------------------------------------
# Make file for UMOA calculation
#--------------------------------------------------
COMPILER=GNU
TARGET=libLinAlg
#--------------------------------------------------
# for compile
#--------------------------------------------------
FDEP=makedepf90
FC=gfortran -ff2c -fPIC # -ff2c option is needed in complex dot product
LIBS=-llapack -lblas
# NOTE -ff2c option is needed in complex dot product

OMP=-fopenmp
FFLAGS=-O3 $(REAL4)
FDFLAGS=-fbounds-check -Wall -fbacktrace -O -Wuninitialized -Ddebug


#--------------------------------------------------
# Source Files
#--------------------------------------------------
SRCDIR = .
DEPDIR = .
SRCC    = $(wildcard $(SRCDIR)/*.c)
SRCF77 += $(wildcard $(SRCDIR)/*.f)
SRCF90 += $(wildcard $(SRCDIR)/*.f90)
SRCF95 += $(wildcard $(SRCDIR)/*.F90)
DEPC = $(SRCC:$(SRCDIR/%.c=$(DEPDIR)/%.d)
SRCS = $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95)

OBJDIR = .
OBJC    = $(SRCC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
OBJF77 += $(SRCF77:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
OBJF90 += $(SRCF90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJF95 += $(SRCF95:$(SRCDIR)/%.F90=$(OBJDIR)/%.o)
OBJS = $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95)

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: $(TARGET)
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -shared -o $(TARGET).so $^ $(LIBS)
	cp $(TARGET).so ../lib/
	cp *.mod ../include/
	@echo '**********************************************************************'
	@echo '* Make sure libLinAlg.so is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
	@echo '**********************************************************************'


$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

dep:
	$(FDEP) $(SRCS) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

clean:
	rm -f $(TARGET).so
	rm -f *.mod *.o

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)


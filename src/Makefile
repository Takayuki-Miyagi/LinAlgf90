#--------------------------------------------------
# Make file for UMOA calculation
#--------------------------------------------------
COMPILER=GNU
TARGET=LinAlg
#--------------------------------------------------
# for compile
#--------------------------------------------------
FDEP=makedepf90
FC=gfortran -ff2c# -ff2c option is needed in complex dot product
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
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $(TARGET).exe $^ $(LIBS)

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
	rm -f $(TARGET).exe
	rm -f *.mod *.o

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)


#--------------------------------------------------
# Make file for UMOA calculation
#--------------------------------------------------
COMPILER=GNU
TARGET=libLinAlg
INSTLDIR=~/LinAlgLib
#--------------------------------------------------
# for compile
#--------------------------------------------------
FDEP=makedepf90
FC=gfortran -ff2c -fPIC
# NOTE -ff2c option is needed in complex dot product
LIBS=-llapack -lblas

OMP=-fopenmp
FFLAGS=-O3
FDFLAGS=-fbounds-check -Wall -fbacktrace -O -Wuninitialized -Ddebug # For debug

#--------------------------------------------------
# Source Files
#--------------------------------------------------
SRCDIR = src
DEPDIR = .
SRCF90 += $(wildcard $(SRCDIR)/*.f90)
SRCF95 += $(wildcard $(SRCDIR)/*.F90)
DEPC = $(SRCC:$(SRCDIR/%.c=$(DEPDIR)/%.d))
SRCS = $(SRCF90) $(SRCF95)

MODDIR = mod
MODF90 += $(SRCF90:$(SRCDIR)/%.f90=$(MODDIR)/%.mod)
MODF95 += $(SRCF95:$(SRCDIR)/%.F90=$(MODDIR)/%.mod)
MODSUP = $(MODF90) $(MODF95)
MODS = $(shell echo $(MODSUP) | tr A-Z a-z)
#$(info $(MODS))

OBJDIR = obj
OBJF90 += $(SRCF90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJF95 += $(SRCF95:$(SRCDIR)/%.F90=$(OBJDIR)/%.o)
OBJS = $(OBJF90) $(OBJF95)

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -shared -o $(TARGET).so $^ $(LIBS)
	if test -d $(INSTLDIR); then \
		: ; \
	else \
		mkdir -p $(INSTLDIR); \
	fi
	ln -sf $(PWD)/$(TARGET).so $(INSTLDIR)/$(TARGET).so
	if test -d $(INSTLDIR)/$(MODDIR); then \
		: ; \
	else \
		mkdir -p $(INSTLDIR)/$(MODDIR); \
	fi
	@for x in $(MODS); do \
		echo linking $(PWD)/$$x '=>' $(INSTLDIR)/$$x; \
		ln -sf $(PWD)/$$x $(INSTLDIR)/$$x; \
	done
	@echo '********************************************************************************'
	@echo '* Make sure $(INSTRIDIR) is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
	@echo '********************************************************************************'

$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(MODDIR)/%.mod:$(SRCDIR)/%.f90 $(OBJDIR)/%.o
	@:

dep:
	$(FDEP) $(SRCS) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi


clean:
	rm -f $(TARGET).so
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJDIR)/*.o

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)


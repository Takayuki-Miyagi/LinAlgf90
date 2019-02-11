#--------------------------------------------------
# Make file for UMOA calculation
#--------------------------------------------------
COMPILER=GNU
TARGET=test_libLinAlg
INSTLDIR=./
#--------------------------------------------------
# for compile
#--------------------------------------------------
FDEP=makedepf90
FC=gfortran
# NOTE -ff2c option is needed in complex dot product
LIBS=-llapack -lblas

OMP=-fopenmp
FFLAGS=-O3
FDFLAGS=-fbounds-check -Wall -fbacktrace -O -Wuninitialized -std=f2008 -Ddebug # For debug

#--------------------------------------------------
# Source Files
#--------------------------------------------------

MAINDIR = test
#SRCMAIN += $(MAINDIR)/test.f90
#SRCMAIN += $(MAINDIR)/test_SVec.f90
#SRCMAIN += $(MAINDIR)/test_DVec.f90
#SRCMAIN += $(MAINDIR)/test_CVec.f90
#SRCMAIN += $(MAINDIR)/test_SMat.f90
#SRCMAIN += $(MAINDIR)/test_DMat.f90
SRCMAIN += $(MAINDIR)/test_CMat.f90
#SRCMAIN += $(MAINDIR)/init_time.f90
#SRCMAIN += $(MAINDIR)/test_type_conv.f90
#SRCMAIN += $(MAINDIR)/test_EigenSolSymD.f90
#SRCMAIN += $(MAINDIR)/test_EigenSolHermite.f90


SRCDIR = src
DEPDIR = .
SRCF90 += $(wildcard $(SRCDIR)/*.f90)
SRCF95 += $(wildcard $(SRCDIR)/*.F90)
SRCS = $(SRCF90) $(SRCF95) $(MAINSRC)

MODDIR = src

OBJDIR = obj
OBJF90 += $(SRCF90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJF95 += $(SRCF95:$(SRCDIR)/%.F90=$(OBJDIR)/%.o)
OBJMAIN += $(SRCMAIN:$(MAINDIR)/%.f90=$(OBJDIR)/%.o)
OBJS = $(OBJF90) $(OBJF95) $(OBJMAIN)

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $(TARGET).exe $^ $(LIBS)
	if test -d $(INSTLDIR); then \
		: ; \
	else \
		mkdir -p $(INSTLDIR); \
	fi
#	ln -sf $(PWD)/$(TARGET).so $(INSTLDIR)/$(TARGET).so
#	if test -d $(INSTLDIR)/$(MODDIR); then \
#		: ; \
	else \
		mkdir -p $(INSTLDIR)/$(MODDIR); \
	fi
#	@for x in $(MODS); do \
#		echo linking $(PWD)/$$x '=>' $(INSTLDIR)/$$x; \
#		ln -sf $(PWD)/$$x $(INSTLDIR)/$$x; \
#	done
#	@echo '********************************************************************************'
#	@echo '* Make sure $(INSTRIDIR) is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
#	@echo '********************************************************************************'

$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -ff2c $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCMAIN)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<

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


# ...
# ...   Title:  Makefile for Transport
# ...
#
#

MODULE := $(TOPD)$(O_DIR)libTransport.a

include $(TOPD)$(TARG_PREFIX).make.cf.$(TARG_NAME)

CURRENT_DIR := ./Transport                     # current directory from $(TOPD)
###############################################################################
# Files

# Source files
wFSRC     :=
FSRC      :=
cSRC      :=
CXXSRC    := src/eqn.cpp src/io.cpp src/LDG.cpp src/limiter.cpp \
             src/problem.cpp src/quad.cpp src/transport.cpp src/utilities.cpp

override CXXFLAGS := $(CXXFLAGS) -Iinclude

include $(TOPD)make.mk


$(MODULE)(%.o) : src/%.cpp
	$(COMPILE.cpp) $<
	$(AR) $(ARFLAGS) $@ $%
	$(RM) $*.o


#------------------------------------------------------------------------
# Dependency List:

ifeq (.depend,$(wildcard .depend))
  include .depend
endif

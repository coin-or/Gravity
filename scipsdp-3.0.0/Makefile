#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
#/* semidefinite programs based on SCIP.                                      */
#/*                                                                           */
#/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
#/*                         EDOM, FAU Erlangen-Nürnberg                       */
#/*               2014-2017 Discrete Optimization, TU Darmstadt               */
#/*                                                                           */
#/*                                                                           */
#/* This program is free software; you can redistribute it and/or             */
#/* modify it under the terms of the GNU Lesser General Public License        */
#/* as published by the Free Software Foundation; either version 3            */
#/* of the License, or (at your option) any later version.                    */
#/*                                                                           */
#/* This program is distributed in the hope that it will be useful,           */
#/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
#/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
#/* GNU Lesser General Public License for more details.                       */
#/*                                                                           */
#/* You should have received a copy of the GNU Lesser General Public License  */
#/* along with this program; if not, write to the Free Software               */
#/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
#/*                                                                           */
#/*                                                                           */
#/* Based on SCIP - Solving Constraint Integer Programs                       */
#/* Copyright (C) 2002-2017 Zuse Institute Berlin                             */
#/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
#/* see file COPYING in the SCIP distribution.                                */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#@file    Makefile
#@brief   Makefile for C++ SDP-Interface for SCIP
#@author  Sonja Mars
#@author  Lars Schewe
#@author  Marc Pfetsch
#@author  Tristan Gally
#@author  Ambros Gleixner

#-----------------------------------------------------------------------------
# own variables
#-----------------------------------------------------------------------------

SCIPSDPVERSION	=	3.0.0
SCIPSDPGITHASH	=
SDPS		=	none

GCCWARN		+= 	-Wextra

# parameters
TIME     	=  	3600
NODES           =       2100000000
MEM		=	6144
THREADS         =       1
PERMUTE         =       0
DISPFREQ	=	10000
FEASTOL		=	default
TEST		=	short
SETTINGS        =       default
CONTINUE	=	false
LOCK		=	false
VALGRIND	=	false
CLIENTTMPDIR    =       /tmp
OPTCOMMAND	=	optimize
SOFTLINKS	=
MAKESOFTLINKS	=	true
OPENBLAS	=	true


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

# save directory to be able to locate library files
ifeq ($(OSTYPE),mingw)
SCIPSDPDIR	=	./
else
SCIPSDPDIR	=	$(realpath .)
endif
SCIPSDPLIBDIR	=	lib

SCIPDIR		= 	$(SCIPSDPDIR)/lib/scip

# check whether SCIPDIR exists
ifeq ("$(wildcard $(SCIPDIR))","")
$(error Please add a soft-link to the SCIP directory as $(SCIPSDPDIR)/lib/scip)
endif

# load SCIP project makefile
include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# settings for SDP solver
#-----------------------------------------------------------------------------

LDFLAGS 	+= 	-lobjscip

SDPIOPTIONS	=
SDPIINC		=
SDPILIB		=
SDPICSRC	=
SDPICCSRC	=
SDPICCOBJ	=


#-----------------------------------------------------------------------------
# DSDP solver
SDPIOPTIONS	+=	dsdp
ifeq ($(SDPS),dsdp)
SDPILIB		= 	$(SCIPSDPLIBDIR)/static/libdsdp.$(STATICLIBEXT) -llapack -lblas
SDPIINC		= 	-I$(SCIPSDPLIBDIR)/include/dsdpinc
SDPICSRC 	= 	src/sdpi/sdpisolver_dsdp.c \
			src/sdpi/lapack_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_dsdp.o \
			$(OBJDIR)/sdpi/lapack_dsdp.o
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/dsdpinc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libdsdp.$(STATICLIBEXT)
SDPIINSTMSG	=	" -> \"dsdpinc\" is the path to the DSDP \"include\" directory, e.g., \"<DSDP-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libdsdp.*\" is the path to the DSDP library, e.g., \"<DSDP-path>/lib/libdsdp.$(STATICLIBEXT)\""
endif

#-----------------------------------------------------------------------------
# SDPA solver
SDPIOPTIONS	+=	sdpa
ifeq ($(SDPS),sdpa)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/sdpainc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libsdpa.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libsdpa.$(SHAREDLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/mumpsinc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libdmumps.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libmumps_common.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libpord.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libmpiseq.$(STATICLIBEXT)
ifeq ($(OPENBLAS),true)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libopenblas.$(SHAREDLIBEXT).0
endif
SDPIINSTMSG	=	" -> \"sdpainc\" is the path to the SDPA \"include\" directory, e.g., \"<SDPA-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libsdpa.*\" is the path to the SDPA library, e.g., \"<SDPA-path>/lib/libsdpa.a\".\n"
SDPIINSTMSG	=	" -> \"mumpsinc\" is the path to the mumps \"include\" directory, e.g., \"<SDPA-path>/mumps/include\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps.*\" is the path to the dmumps library, e.g., \"<SDPA-path>/mumps/build/lib/libdmumps.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps_common.*\" is the path to the mumps_common library, e.g., \"<SDPA-path>/mumps/build/lib/libmumps_common.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps_common.*\" is the path to the pord library, e.g., \"<SDPA-path>/mumps/build/lib/libpord.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps.*\" is the path to the mpiseq library, e.g., \"<SDPA-path>/mumps/build/libseq/libmpiseq.$(STATICLIBEXT)\".\n"
ifeq ($(OPENBLAS),true)
SDPIINSTMSG	+=	" -> \"libopenblas.$(SHAREDLIBEXT).0\" is the openblas library.\n"
SDPILIB		=      	-L$(SCIPSDPLIBDIR)/$(LIBEXTTYPE) -lsdpa $(SCIPSDPLIBDIR)/static/libdmumps.$(STATICLIBEXT) $(SCIPSDPLIBDIR)/static/libmumps_common.$(STATICLIBEXT) \
			$(SCIPSDPLIBDIR)/static/libpord.$(STATICLIBEXT)	$(SCIPSDPLIBDIR)/static/libmpiseq.$(STATICLIBEXT) $(SCIPSDPLIBDIR)/shared/libopenblas.$(SHAREDLIBEXT).0 \
			-Wl,-rpath,$(SCIPSDPDIR)/$(SCIPSDPLIBDIR)/static -Wl,-rpath,$(SCIPSDPDIR)/$(SCIPSDPLIBDIR)/shared -lgfortran -L/lib/x86_64-linux-gnu -lpthread -lgomp
else
SDPILIB		=      	-L$(SCIPSDPLIBDIR)/$(LIBEXTTYPE) -lsdpa $(SCIPSDPLIBDIR)/static/libdmumps.$(STATICLIBEXT) $(SCIPSDPLIBDIR)/static/libmumps_common.$(STATICLIBEXT) \
			$(SCIPSDPLIBDIR)/static/libpord.$(STATICLIBEXT)	$(SCIPSDPLIBDIR)/static/libmpiseq.$(STATICLIBEXT) -lgfortran -llapack -lblas
endif
SDPIINC	=  -I$(SCIPSDPLIBDIR)/include/sdpainc
SDPIINC	+= -I$(SCIPSDPLIBDIR)/include/mumpsinc
SDPICCSRC 	= 	src/sdpi/sdpisolver_sdpa.cpp
SDPICSRC	=	src/sdpi/lapack_sdpa.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_sdpa.o $(OBJDIR)/sdpi/lapack_sdpa.o
endif

ifneq ($(SDPS),sdpa)
DISABLEOMP=1
endif

#-----------------------------------------------------------------------------
# MOSEK solver
SDPIOPTIONS	+=	msk
ifeq ($(SDPS),msk)
# decide on 32 or 64 bit
BITEXT     =  $(word 2, $(subst _, ,$(ARCH)))
SDPIINC		= 	-I$(SCIPSDPLIBDIR)/include/mosekh
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/mosekh
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libmosek$(BITEXT).$(SHAREDLIBEXT)
SDPIINSTMSG	=	"  -> \"mosekh\" is the path to the MOSEK \"h\" directory, e.g., \"<MOSEK-path>/8/tools/platform/linux64x86/h\".\n"
SDPIINSTMSG	+=	" -> \"libmosek$(BITEXT).*\" is the path to the MOSEK library, e.g., \"<MOSEK-path>/8/tools/platform/linux64x86/bin/libmosek$(BITEXT).$(SHAREDLIBEXT)\".\n"
SDPILIB		= 	-m$(BITEXT) $(SCIPSDPLIBDIR)/shared/libmosek$(BITEXT).$(SHAREDLIBEXT) -Wl,-rpath=$(dir $(realpath $(SCIPSDPDIR)/$(SCIPSDPLIBDIR)/shared/libmosek$(BITEXT).$(SHAREDLIBEXT))) -llapack -lblas -pthread -lc -lm
SDPICSRC 	= 	src/sdpi/sdpisolver_mosek.c src/sdpi/lapack_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_mosek.o $(OBJDIR)/sdpi/lapack_dsdp.o
endif

#-----------------------------------------------------------------------------
# no solver
SDPIOPTIONS	+=	none
ifeq ($(SDPS),none)
SDPILIB		= 	-L$(SCIPSDPLIBDIR) -llapack -lblas
SDPICSRC 	= 	src/sdpi/sdpisolver_none.c src/sdpi/lapack_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_none.o $(OBJDIR)/sdpi/lapack_dsdp.o
SETTINGS	= 	lp_approx
endif

# include install/uninstall targets
-include make/make.install

LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT).$(GAMS)

#-----------------------------------------------------------------------------

SDPOBJSUBDIRS	=	$(OBJDIR)/scipsdp \
			$(OBJDIR)/sdpi

#-----------------------------------------------------------------------------
# OMPSETTINGS (used to set number of threads for Openblas)
#-----------------------------------------------------------------------------

ifeq ($(OMP),false)
DISABLEOMP=1
endif

ifeq ($(OMP),true)
DISABLEOMP=0
endif

OMPFLAGS =
ifneq ($(DISABLEOMP),1)
OMPFLAGS += -DOMP
endif

#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME	=	scipsdp
MAINCOBJ	=	scipsdp/SdpVarmapper.o \
			scipsdp/SdpVarfixer.o \
			scipsdp/cons_sdp.o \
			scipsdp/cons_savedsdpsettings.o \
			scipsdp/relax_sdp.o \
			scipsdp/disp_sdpiterations.o \
			scipsdp/disp_sdpavgiterations.o \
			scipsdp/disp_sdpfastsettings.o \
			scipsdp/disp_sdppenalty.o \
			scipsdp/disp_sdpunsolved.o \
			scipsdp/prop_sdpredcost.o \
			scipsdp/branch_sdpmostfrac.o \
			scipsdp/branch_sdpmostinf.o \
			scipsdp/branch_sdpobjective.o \
			scipsdp/branch_sdpinfobjective.o \
			scipsdp/heur_sdpfracdiving.o \
			scipsdp/heur_sdprand.o \
			scipsdp/reader_cbf.o \
			scipsdp/prop_sdpobbt.o \
			sdpi/sdpi.o \
			sdpi/sdpsolchecker.o \
			scipsdpgithash.o

MAINCCOBJ 	=	scipsdp/main.o \
			scipsdp/objreader_sdpa.o \
			scipsdp/ScipStreamBuffer.o


MAINCSRC	=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))
MAINCCSRC 	=	$(addprefix $(SRCDIR)/,$(MAINCCOBJ:.o=.cpp))
MAINDEP 	=	$(SRCDIR)/depend.cppmain.$(OPT)

SCIPSDPGITHASHFILE	= 	$(SRCDIR)/scipsdpgithash.c

# @todo possibly add LPS
MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(SDPS)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCOBJ))
MAINCCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCCOBJ))

ALLSRC		=	$(MAINCSRC) $(MAINCCSRC) $(SDPICSRC) $(SDPICCSRC)
LINKSMARKERFILE =	$(SCIPSDPLIBDIR)/linkscreated.$(SDPS).$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LASTSETTINGS 	=	$(OBJDIR)/make.lastsettings


#-----------------------------------------------------------------------------
# rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(SDPIOBJ) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: checkdefines
checkdefines:
ifeq ($(SDPIOBJ),)
		$(error invalid SDP solver selected: SDPS=$(SDPIS). Possible options are: $(SDPIOPTIONS))
endif

.PHONY: preprocess
preprocess:     checkdefines
		@$(SHELL) -ec 'if test ! -e $(LINKSMARKERFILE) ; \
			then \
				echo "-> generating necessary links" ; \
				$(MAKE) -j1 $(LINKSMARKERFILE) ; \
			fi'
		@$(MAKE) touchexternal

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*/*.c src/*/*.cpp src/*/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h;

# include target to detect the current git hash
-include make/make.detectgithash

# this empty target is needed for the SCIP release versions
githash::      # do not remove the double-colon

.PHONY: lint
lint:		$(MAINCSRC) $(MAINCCSRC) $(SDPICSRC) $(SDPICCSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -i lint co-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK): $(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@mkdir -p $(OBJDIR);

$(SDPOBJSUBDIRS):	| $(OBJDIR)
		@-mkdir -p $(SDPOBJSUBDIRS);

$(SCIPSDPLIBDIR):
		@-mkdir -p $(SCIPSDPLIBDIR)
		
$(SCIPSDPLIBDIR)/include: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/include
		
$(SCIPSDPLIBDIR)/static: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/static
		
$(SCIPSDPLIBDIR)/shared: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/shared

$(BINDIR):
		-@test -d $(BINDIR) || { \
		echo "-> Creating $(BINDIR) directory"; \
		mkdir -p $(BINDIR); }

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rm -f $(OBJDIR)/scipsdp/*.o
		@-rm -f $(OBJDIR)/sdpi/*.o
		@-rm -f $(OBJDIR)/*.o
		@-rmdir $(OBJDIR)/scipsdp
	 	@-rmdir $(OBJDIR)/sdpi
		@-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

#-----------------------------------------------------------------------------
-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(SCIPSDPLIBDIR) $(OBJDIR)
ifneq ($(USRFLAGS),$(LAST_USRFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USROFLAGS),$(LAST_USROFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCFLAGS),$(LAST_USRCFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRLDFLAGS),$(LAST_USRLDFLAGS))
		@-touch -c $(ALLSRC)
endif
ifneq ($(USRARFLAGS),$(LAST_USRARFLAGS))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBLKMEM),$(LAST_NOBLKMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBUFMEM),$(LAST_NOBUFMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBLKBUFMEM),$(LAST_NOBLKBUFMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(SDPS),$(LAST_SDPS))
		@-touch -c $(SRCDIR)/scipsdp/cons_sdp.c
endif
ifneq ($(OMP),$(LAST_OMP))
		@-touch -c $(SRCDIR)/scipsdp/cons_sdp.c
endif
ifneq ($(SCIPSDPGITHASH),$(LAST_SCIPSDPGITHASH))
		@$(MAKE) githash
endif
		@$(SHELL) -ec 'if test ! -e $(SCIPGITHASHFILE) ; \
			then \
				echo "-> generating $(SCIPGITHASHFILE)" ; \
				$(MAKE) githash ; \
			fi'
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_PARASCIP=$(PARASCIP)" >> $(LASTSETTINGS)
		@echo "LAST_USRFLAGS=$(USRFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USROFLAGS=$(USROFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRDFLAGS=$(USRDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKMEM=$(NOBLKMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBUFMEM=$(NOBUFMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKBUFMEM=$(NOBLKBUFMEM)" >> $(LASTSETTINGS)
		@echo "LAST_SDPS=$(SDPS)" >> $(LASTSETTINGS)
		@echo "LAST_OMP=$(OMP)" >> $(LASTSETTINGS)
		@echo "LAST_SCIPGITHASH=$(SCIPSDPGITHASH)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		| $(SCIPSDPLIBDIR) $(SCIPSDPLIBDIR)/include $(SCIPSDPLIBDIR)/static $(SCIPSDPLIBDIR)/shared echosoftlinks $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: SDPS=$(SDPS) LPS=$(LPS) SUFFIX=$(LINKLIBSUFFIX) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP)"
		@echo
		@echo "* SCIPSDP needs some softlinks to external programs, in particular, SDP-solvers."
		@echo "* Please insert the paths to the corresponding directories/libraries below."
		@echo "* The links will be installed in the 'lib' directory."
		@echo "* For more information and if you experience problems see the INSTALL file."
		@echo
		@echo -e $(SDPIINSTMSG)

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				BASENAMEA=`basename $@ .$(STATICLIBEXT)` ; \
				BASENAMESO=`basename $@ .$(SHAREDLIBEXT)` ; \
				echo ; \
				echo "- preparing missing soft-link \"$@\":" ; \
				if test -e $$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				if test -e $$DIRNAME/$$BASENAMESO.$(STATICLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMESO.$(STATICLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(SCIPSDPDIR) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

#-----------------------------------------------------------------------------
.PHONY: test
test:
		@-(cd check && ln -fs $(SCIPDIR)/check/check.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/allcmpres.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs $(SCIPDIR)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPDIR)/check/getlastprob.awk);
		@-(cd check && ln -fs $(SCIPDIR)/check/configuration_set.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/configuration_tmpfile_setup_scip.sh configuration_tmpfile_setup_$(MAINNAME).sh);
		@-(cd check && ln -fs $(SCIPDIR)/check/run.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) \
			$(CONTINUE) $(LOCK) $(SCIPSDPVERSION) $(SDPS) $(DEBUGTOOL) $(CLIENTTMPDIR) $(REOPT) $(OPTCOMMAND) $(SETCUTOFF) $(MAXJOBS) $(VISUALIZE) $(PERMUTE) $(SEEDS);

#-----------------------------------------------------------------------------

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(SDPIINC) $(DFLAGS) $(MAINCCSRC) $(SDPICCSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.cpp|g'\'' \
		>$(MAINDEP)'
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(SDPIINC) $(DFLAGS) $(MAINCSRC) $(SDPICSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.c|g'\'' \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/sdpi/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/sdpi/\2.c|g'\'' \
		>>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	preprocess $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(SDPIOBJ) | $(SDPOBJSUBDIRS) $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(LINKCCSCIPALL) $(SDPIOBJ) $(SDPILIB) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CFLAGS) $(OMPFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CXXFLAGS) $(OMPFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------

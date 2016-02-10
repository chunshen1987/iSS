# ===========================================================================
#  Makefile iSS                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 -Wall $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	iSS.e
endif

SRC		=	main.cpp arsenal.cpp ParameterReader.cpp \
			RandomVariable1DArray.cpp RandomVariable2DArray.cpp \
			RandomVariable.cpp NBD.cpp TableFunction.cpp Table.cpp \
			emissionfunction.cpp readindata.cpp Poisson.cpp

INC		= 	main.h arsenal.h ParameterReader.h \
			RandomVariable1DArray.h RandomVariable2DArray.h \
			RandomVariable.h NBD.h TableFunction.h Table.h \
			emissionfunction.h readindata.h Stopwatch.h \
                  Poisson.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

./NBD.cpp:  arsenal.h RandomVariable.h Table.h TableFunction.h
./Poisson.cpp:  arsenal.h RandomVariable.h Table.h TableFunction.h
./ParameterReader.cpp:  arsenal.h
./RandomVariable.cpp:  arsenal.h TableFunction.h
./RandomVariable1DArray.cpp:  arsenal.h
./RandomVariable2DArray.cpp:  arsenal.h
./RandomVariableNDArray.cpp:  arsenal.h
./Table.cpp:  arsenal.h
./TableFunction.cpp:  arsenal.h
./emissionfunction.cpp: main.h readindata.h RandomVariable1DArray.h RandomVariable2DArray.h RandomVariableNDArray.h NBD.h ParameterReader.h arsenal.h Stopwatch.h Table.h Poisson.h 
./main.cpp: ParameterReader.h Table.h readindata.h emissionfunction.h arsenal.h
./readindata.cpp: main.h arsenal.h Table.h

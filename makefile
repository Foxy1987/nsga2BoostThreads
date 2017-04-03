SHELL = /bin/sh


CC = g++
ROOT=/afs/cad.njit.edu/research/dms/horacio/1
LIBS = -lboost_thread -lboost_system libsimulate.a libzf.a
CPPFLAGS = -g -O2 -Wall

VPATH=%.h ./include
VPATH=%.o ./obj

OBJDIR = ./obj
INCLUDE_DIR = ./include 

INCLUDES  := $(addprefix -I,$(INCLUDE_DIR)) 


objects = $(addprefix $(OBJDIR)/, main.o individual.o costfunc.o population.o random.o worker.o)

MY_APPS = test

$(MY_APPS) : $(objects)
	${CC} -o app1 ${objects} ${CPPFLAGS} ${LIBS} -lm

$(OBJDIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) ${INCLUDES} $< -o $@


$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	rm -f ${MY_APPS}
	rm -f ${objects}

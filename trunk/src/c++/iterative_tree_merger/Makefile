
# ------------ Programs --------------

CC		= gcc
CXX		= g++
CD		= cd
AR		= ar rc
RM		= rm -f
ECHO		= @echo -e
ETAGS		= etags
PRINT		= a2ps
MAKEDEPEND 	= makedepend

# ----------- Object Rules ----------



SOURCES		= 

PROG_SOURCES	= FNJ_collapse_method.cpp iterative_merge.cpp


# -------- Derieved Parameters ------

C_SOURCES	= $(filter %.c, $(SOURCES))
CPP_SOURCES	= $(filter %.cpp, $(SOURCES))
C_PROGS		= $(patsubst %.c, %, $(filter %.c, $(PROG_SOURCES)))
CPP_PROGS	= $(patsubst %.cpp, %, $(filter %.cpp, $(PROG_SOURCES)))



C_OBJECTS	= $(C_SOURCES:%.c=%.o)
CPP_OBJECTS	= $(CPP_SOURCES:%.cpp=%.o)
OBJECTS		= $(C_OBJECTS) $(CPP_OBJECTS)



# -------- Compiler Flags ------------------

INCS		+= -I/usr/include -I. -I../ -I../distance_methods -L../
LIBS		+= -lm

CFLAGS		+= -Wall  $(DEFINES) $(INCS) 
#CFLAGS		+= -ggdb -g
CFLAGS		+= -O3 -DNDEBUG
#CFLAGS		+= -pg

# -------- Targets -------------------------


all		: cleanprogs $(C_OBJECTS) $(CPP_OBJECTS) $(C_PROGS) $(CPP_PROGS)

$(C_OBJECTS): %.o: %.c
		$(CC) -c $(CFLAGS) $< -o $@
$(CPP_OBJECTS): %.o: %.cpp
		$(CXX) -c $(CFLAGS) $< -o $@

$(C_PROGS):  %: %.c $(patsubst %, %, $(OBJECTS))
		$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(CPP_PROGS):  %: %.cpp $(patsubst %, %, $(OBJECTS)) ../compbiolib.ar
		$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)


clean		: cleanprogs
		$(RM) $(patsubst %, %, $(OBJECTS))

cleanprogs	:
		$(RM) $(patsubst %, %, $(C_PROGS) $(CPP_PROGS) core)
		$(RM) $(patsubst %, %.exe, $(C_PROGS) $(CPP_PROGS)) #under Windows


listvars	:
		$(ECHO) C_SOURCES: $(C_SOURCES)
		$(ECHO) C_OBJECTS: $(C_OBJECTS)
		$(ECHO) CPP_SOURCES: $(CPP_SOURCES)
		$(ECHO) CPP_OBJECTS: $(CPP_OBJECTS)
		$(ECHO) C_PROGS: $(C_PROGS)
		$(ECHO) CPP_PROGS: $(CPP_PROGS)

depend		:
		$(MAKEDEPEND) -- $(CFLAGS) -- $(SOURCES)

etags		:
		$(ETAGS) $(SOURCES) $(PROG_SOURCES)
print		:
		$(PRINT) $(SOURCES) $(PROG_SOURCES)

# -------- Dependencies --------------------


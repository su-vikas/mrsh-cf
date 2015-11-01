CC=g++

LIB = $(wildcard lib/cuckoofilter/src/*.h)
PROJECT_SRC = ./src/main.cpp  ./src/fnv.c ./lib/cuckoofilter/src/hashutil.cpp #./lib/cuckoofilter/src/printutil.cpp
DBGF = -g -ggdb -O0

CFLAGS += -std=c++11  -I. -I./include -I/usr/include/boost -I/usr/include -I/usr/local/boost -I./src/ -I/usr/local/include
#-Wall

OPT = -O3  # O3 gives same performance as O2 but O3 might become unpredictable. 


#LDFLAGS += -Wall -lpthread -lssl -lcrypto -lm -lboost_filesystem -lboost_system 
LDFLAGS += -lpthread -lssl -lcrypto -lm -Wl,-rpath=/usr/local/lib -lboost_filesystem -lboost_system -lboost_program_options

# HEADERS 

NAME=mrshcf.exe

all:debug

debug: ${PROJECT_SRC} ${LIB} Makefile
	${CC} $(DBGF) $(CFLAGS)  -o ${NAME} ${PROJECT_SRC} $(LIB) $(LDFLAGS)

mrsh: ${PROJECT_SRC} ${PROJECT_HDR}
	${CC} $(OPT) $(CFLAGS)  -o ${NAME} ${PROJECT_SRC} $(LIB) $(LDFLAGS)
#pg for profiler, gprof. 

prof: ${PROJECT_SRC} ${PROJECT_HDR}
	${CC} $(OPT) $(CFLAGS) -pg  -o ${NAME} ${PROJECT_SRC} $(LIB) $(LDFLAGS)

clean :  
	rm -f mrsh_cuckoo.exe *.o 

#%.o: %.cpp ${HEADERS} Makefile
#	$(CC) $(CFLAGS) $< -o $@

#for DT_DIR feature to work, need to have the _BSD_SOURCE  feature test macro defined. THese are not standard, and GCC does not define the macro when compiling for C99
# -lm: -l means link a library and -m means a math library. Without this option 



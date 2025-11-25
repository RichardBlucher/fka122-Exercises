CC = gcc
LOCAL=${HOME}/.local
LIBDIR=lib

CFLAGS = \
	 -g \
	 -pedantic \
	 -Werror \
	 -Wall \
	 -Iinclude \
	 -Wl,-rpath,${LOCAL} \
	 -I${LOCAL}/include \
	 -Wl,-rpath,${LIBDIR}\
	 -L${LIBDIR} \
	 -L${LOCAL}/lib \
	 -Isrc/nepc/



CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas \
      -lnep

OBJ = \
      obj/tools.o \
      obj/lattice.o \
      obj/potential.o \
      obj/run.o

OBJ_LIB = \
	obj_nep/nep.o \
	obj_nep/interface.o

CXXFLAGS = \
	   -DUSE_TABLE_FOR_RADIAL_FUNCTIONS \
	   -O3 \
	   -fPIC

MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

_H1: obj lib H1

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@

H1: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ $(LIB) -o $@

obj_nep/%.o: src/nepc/%.cpp
	$(CXX) -MMD -c $(CXXFLAGS) $< -o $@

lib: $(OBJ_LIB)
	$(CXX) -shared $^ -o ${LIBDIR}/libnep.so

obj:
	mkdir -p obj
	mkdir -p obj_nep
	mkdir -p ${LIBDIR}

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f H1 run-test
	rm -rf obj obj_test obj_nep ${LIBDIR}

.PHONY: clean

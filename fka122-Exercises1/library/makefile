CC = gcc

LOCAL = ${HOME}/.local
CFLAGS = \
	 -pedantic \
	 -Wall \
	 -Werror \
	 -g \
	 -Iinclude \
	 -MMD \
	 -Wl,-rpath,${LOCAL}/lib \
	 -L${LOCAL}/lib \
	 -I${LOCAL}/include


CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

program: obj _program

_program: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o program $(LIB)

obj/%.o: src/%.c
	$(CC) -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test

.PHONY: clean

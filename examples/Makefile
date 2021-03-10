all: test_complexe2 test_complexe test_dot test_dot_dyn test_swap test_copy test_amax test_amin test_axpy

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -fopenmp -L$(LIB_DIR) $(LIBST)  
OPTIONS_LINK_DYN  = -fopenmp -L$(LIB_DIR) $(LIBDYN)  

test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

test_dot_dyn: test_dot.o flop.o
	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN)

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_dot.c

test_complexe.o: test_complexe.c  $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_complexe: test_complexe.o flop.o
	gcc -o test_complexe test_complexe.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe2.o: test_complexe2.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe2.c

test_complexe2: test_complexe2.o flop.o
	gcc -o test_complexe2 test_complexe2.o flop.o $(OPTIONS_LINK_STATIC)

test_swap.o: test_swap.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_swap.c

test_swap: test_swap.o $(LIB_DIR)/libmnblas.a
	gcc -o test_swap test_swap.o  $(OPTIONS_LINK_STATIC)

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_copy: test_copy.o $(LIB_DIR)/libmnblas.a
	gcc -o test_copy test_copy.o  $(OPTIONS_LINK_STATIC)

test_amax.o: test_amax.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_amax.c

test_amax: test_amax.o $(LIB_DIR)/libmnblas.a
	gcc -o test_amax test_amax.o  $(OPTIONS_LINK_STATIC)

test_amin.o: test_amin.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_amin.c

test_amin: test_amin.o $(LIB_DIR)/libmnblas.a
	gcc -o test_amin test_amin.o  $(OPTIONS_LINK_STATIC)

test_axpy.o: test_axpy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_axpy.c

test_axpy: test_axpy.o $(LIB_DIR)/libmnblas.a
	gcc -o test_axpy test_axpy.o  $(OPTIONS_LINK_STATIC)

clean:
	rm -f *.o test_dot test_dot_dyn test_complexe test_complexe2 *~


OBJECTS = main.o int_three_center.o int_overlap.o int_kinetic.o int_nuclear.o int_eri.o igamma.o

OPTIMIZE = -O3
LIBRARY  = -lm

test:	${OBJECTS}
	cc ${OPTIMIZE} -o test.x ${OBJECTS} ${LIBRARY}

clean:
	rm *.o test.x

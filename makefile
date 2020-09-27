EXECS=transportation_sort
MPICC?=mpicc

all: ${EXECS}

transportation_sort: transportation_sort.c
	${MPICC} -o transportation_sort transportation_sort.c

clean:
	rm -f ${EXECS}
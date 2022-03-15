CC = mpicc
RM = /bin/rm
PROG = mm mm_abft
CFLAGS = -lm -lz

DIR_SRC = ./src
DIR_BIN = ./bin
DIR_MM = ${DIR_SRC}/mm
LIB_ABFT = ${DIR_SRC}/abft.c

all : ${PROG}

mm : ${DIR_BIN}/mm
${DIR_BIN}/mm : ${DIR_MM}/mm.c 
	${CC} -o $@ $< ${CFLAGS}

mm_abft : ${DIR_BIN}/mm_abft
${DIR_BIN}/mm_abft : ${DIR_MM}/mm_abft.c ${LIB_ABFT}
	${CC} -o $@ $^ ${CFLAGS}

clean :
	${RM} -f ${DIR_BIN}/*







GSL_LIBS=-lgsl -lgslcblas
GSL_INC=
FFTW_LIBS=-lfftw3_mpi -lfftw3
FFTW_INC=

GLIBPREFIX=../../GLIB-STATIC
CFLAGS = -I.. -I. -I$(GLIBPREFIX)/include/glib-2.0 -I$(GLIBPREFIX)/lib/glib-2.0/include $(GSL_INC) $(FFTW_INC)
LDFLAGS = -L$(GLIBPREFIX)/lib
CC = mpiicc -g

SOURCES = \
    main.c \
    ../mpiu.c \
    commonblock.c \
    power.c \
    disp.c \
	filter.c

HEADERS =
OBJECTS = $(SOURCES:.c=.o)
.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

ngenic.grid: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lglib-2.0 -lm $(GSL_LIBS) $(FFTW_LIBS)

clean: 
	rm -rf $(OBJECTS) ngenic.grid

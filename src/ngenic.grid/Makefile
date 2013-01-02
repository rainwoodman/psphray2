include Makefile.local

CFLAGS = -I.. -I. $(GLIB_INC) $(GSL_INC) $(FFTW_INC)
LDFLAGS =

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
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(GLIB_LIBS) -lm $(GSL_LIBS) $(FFTW_LIBS)

clean: 
	rm -rf $(OBJECTS) ngenic.grid
Makefile.local: Makefile.local.template
	if ! [ -f $@ ]; then cp $< $@; else echo Makefile.local exsits; fi

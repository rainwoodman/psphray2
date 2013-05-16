include Makefile.local

CFLAGS = -I../src -I. $(GLIB_INC) $(GSL_INC) $(FFTW_INC)
LDFLAGS =

SOURCES = \
    main.c \
    ../src/mpiu.c \
    commonblock.c \
    power.c \
	degrade.c \
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
	@if ! [ -f $@ ]; then cp $< $@; else if ! diff -y $< $@; then echo please edit and save Makefile.local; fi; fi

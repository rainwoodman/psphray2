include Makefile.local

CFLAGS = -I. -I../src $(GLIB_INC) $(GSL_INC) $(FFTW_INC)
LDFLAGS =

SOURCES = \
    ../src/paramfile.c \
    ../src/commonblock.c \
	snapshot.c \
    mpiu.c \
    par.c \
	domain.c 
TESTS = tests/test_par.c

HEADERS =
OBJECTS = $(SOURCES:.c=.o)
.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

attic: $(OBJECTS) main.c
	$(CC) -o $@ main.c $(OBJECTS) $(CFLAGS) $(LDFLAGS) $(GLIB_LIBS) -lm $(GSL_LIBS) $(FFTW_LIBS)

test:  $(TESTS)
	for test in $(TESTS); do \
		( \
		$(CC) $(CFLAGS) -o $${test//.c} $(OBJECTS) $${test} $(LDFLAGS) $(GLIB_LIBS) -lm $(GSL_LIBS) $(FFTW_LIBS); \
		) \
	done;

clean: 
	rm -rf $(OBJECTS) attic

Makefile.local: Makefile.local.template
	@if ! [ -f $@ ]; then cp $< $@; else if ! diff -y $< $@; then echo please edit and save Makefile.local; fi; fi

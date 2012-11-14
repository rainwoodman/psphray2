GLIBPREFIX=../GLIB-STATIC
CFLAGS = -I$(GLIBPREFIX)/include/glib-2.0 -I$(GLIBPREFIX)/lib/glib-2.0/include
LDFLAGS = -L$(GLIBPREFIX)/lib
CC = mpiicc -g

SOURCES = \
    main.c \
    commonblock.c \
    par.c \
    paramfile.c \
    snapshot.c

HEADERS =
OBJECTS = $(SOURCES:.c=.o)
.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

psphray2: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lglib-2.0 -lm

clean: 
	rm -rf $(OBJECTS) psphray2

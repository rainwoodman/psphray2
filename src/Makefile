GLIBPREFIX=../GLIB-STATIC
CFLAGS = -I. -I$(GLIBPREFIX)/include/glib-2.0 -I$(GLIBPREFIX)/lib/glib-2.0/include
LDFLAGS = -L$(GLIBPREFIX)/lib
CC = mpiicc -g

GADGETSOURCES = gadget/gadgetall.c gadget/cooling.c gadget/sfr_eff.c
SOURCES = \
    main.c \
    commonblock.c \
    fckey.c \
    par.c \
    paramfile.c \
    snapshot.c \
	domain.c \
	tree.c \
	$(GADGETSOURCES)

HEADERS =
OBJECTS = $(SOURCES:.c=.o)
.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

psphray2: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lglib-2.0 -lm

clean: 
	rm -rf $(OBJECTS) psphray2

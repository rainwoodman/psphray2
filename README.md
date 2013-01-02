psphray2
========

Parallel SPH Ray Tracing v2

you will need to install a staticly build GLib-2 to the directory GLIB-STATIC.
Note: GLib-2 has a lot of stuff (glib, gthread, gmodule, gobject) and we only 
use the core module (glib).  

To compile the static version of glib into GLIB-STATE, at the toplevel code
directory, do the followings:

wget http://ftp.gnome.org/pub/gnome/sources/glib/2.30/glib-2.30.2.tar.bz2
tar -xjvf glib-2.30.2.tar.bz2
cd glib-2.30.2
./configure --prefix=$PWD/../GLIB-STATIC --enable-static --disable-shared
cd glib
make install

There is a GADGET IC generator at src/ngenic.grid. use make to build it.
See README there.



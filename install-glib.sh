#! bash

wget http://ftp.gnome.org/pub/gnome/sources/glib/2.30/glib-2.30.2.tar.bz2

tar -xjvf glib-2.30.2.tar.bz2
(
  cd glib-2.30.2
  ./configure --prefix=$PWD/../GLIB-STATIC --enable-static --disable-shared --enable-debug=yes LIBFFI_CFLAGS=-DFFI LIBFFI_LIBS=-L/usr/lib64
  (cd glib; make install )
)

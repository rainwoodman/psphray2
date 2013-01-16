#! bash

wget http://ftp.gnome.org/pub/gnome/sources/glib/2.30/glib-2.30.2.tar.bz2

tar -xjvf glib-2.30.2.tar.bz2
(
  cd glib-2.30.2
  ./configure --prefix=$PWD/../GLIB-STATIC --enable-static --disable-shared
  (cd glib; make install )
)

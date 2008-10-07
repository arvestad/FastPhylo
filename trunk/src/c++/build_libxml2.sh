#!/bin/sh
set -e
libxml2version=$1
installDir="$2"
tar xfz "../download/libxml2-${libxml2version}.tar.gz"
cd libxml2-${libxml2version}
./configure  "--prefix=$installDir" "--with-zlib=$installDir/lib/" --without-iconv --without-threads --disable-ipv6 --without-python --without-html --without-ftp --without-http --without-schematron --disable-shared && make libxml2.la && install -D .libs/libxml2.a "$installDir/lib/libxml2.a" && cd include && make install

#!/bin/sh
set -e
libxml2version=$1
installDir="$2"
tar xfz "../download/libxml2-${libxml2version}.tar.gz"
cd libxml2-${libxml2version}
./configure  "--prefix=$installDir" "--with-zlib=$installDir" --without-iconv --without-threads --disable-ipv6 --without-python --without-html --without-ftp --without-http --without-schematron --disable-shared && make && make install

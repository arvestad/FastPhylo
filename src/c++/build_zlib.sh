#!/bin/sh
set -e
zlibversion=$1
installDir="$2"
tar xfz "../download/zlib-${zlibversion}.tar.gz"
cd zlib-${zlibversion}
./configure  "--prefix=$installDir" && make && make install

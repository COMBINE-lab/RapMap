#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

host=
version=

while getopts "v:n:" opt; do
  case $opt in
    n)
      echo "Host is $OPTARG" >&2
      host=$OPTARG
      ;;
    v)
      echo "Version is $OPTARG" >&2
      version=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

echo -e "Preparing binary release\n=====================\n"
echo -e "Version = ${version}"
echo -e "Host = ${host}"

# create the binary directory 
relname=RapMap-${version}_${host}
mkdir -p ${DIR}/../RELEASES
mkdir -p ${DIR}/../RELEASES/${relname}
mkdir -p ${DIR}/../RELEASES/${relname}/bin
mkdir -p ${DIR}/../RELEASES/${relname}/lib

echo -e "Copying over the binary\n"
cp ${DIR}/../bin/rapmap ${DIR}/../RELEASES/${relname}/bin/

# copy other dependencies (shared libraries)
#echo -e "Copying over other shared library dependencies\n"
#bash ${DIR}/../scripts/cpld.bash ${DIR}/../bin/salmon ${DIR}/../RELEASES/${relname}/lib/
#echo -e "Removing dangerous dependencies\n"
#rm ${DIR}/../RELEASES/${relname}/lib/libc.so.6
#rm ${DIR}/../RELEASES/${relname}/lib/ld-linux-x86-64.so.2
#rm ${DIR}/../RELEASES/${relname}/lib/libdl.so.2
#rm ${DIR}/../RELEASES/${relname}/lib/libpthread*.so.*

# now make the tarball
echo -e "Making the tarball\n"
cd ${DIR}/../RELEASES
tar czvf ${relname}.tar.gz ${relname}

echo -e "Done making release!"

#!/bin/bash
set -e
# from http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
no_native_arch=false
cxxflags=""

while [[ $# > 1 ]]
do
    key="$1"

    case $key in
        -b|--branch)
            branch="$2"
            shift # past argument
            ;;
        -v|--version)
            version="$2"
            shift # past argument
            ;;
        -f|--cxxflags)
            cxxflags="$2"
            shift # past argument
            ;;
        --no-native)
            no_native_arch=true
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done

echo "Building rapmap [branch = ${branch}]. Tagging version as ${version}"
if [ "$no_native_arch" = true ] ; then 
    echo "Disabling -march=native"
fi

if [[ -z $cxxflags ]] ; then
    echo "Passed CXXFLAGS ${cxxflags}"
fi

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -x

# Install things we need
yum install -y --quiet wget
wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
rpm -i --quiet epel-release-5-4.noarch.rpm
#yum install -y --quiet git
#yum install -y --quiet xz-devel.x86_64
#yum install -y --quiet bzip2-devel.x86_64
yum install -y --quiet unzip

curl -k -L https://github.com/COMBINE-lab/RapMap/archive/${branch}.zip -o ${branch}.zip
unzip ${branch}.zip
mv RapMap-${branch} RapMap
cd RapMap
mkdir build
cd build

if [ "$no_native_arch" = true ] ; then 
    cmake -DFETCH_BOOST=TRUE -DCMAKE_CXX_FLAGS=${cxxflags} -DNO_NATIVE_ARCH ..
else
    cmake -DFETCH_BOOST=TRUE -DCMAKE_CXX_FLAGS=${cxxflags} ..
fi

make
make install
make test
cd ../scripts
bash make-release.sh -v ${version} -n CentOS5
cd ../RELEASES
cp *.tar.gz /io/

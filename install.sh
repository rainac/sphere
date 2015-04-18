
export PREFIX=$HOME/sw/sph-new

mkdir $PREFIX

sudo apt-get install make g++ libhdf5-dev libvtk5-dev freeglut3-dev octave

wget https://codeforge.lbl.gov/frs/download.php/170/H5Part-1.6.1.tar.gz

tar -xzf H5Part-1.6.1.tar.gz

pushd /usr/local
mkdir hdf5
cd hdf5
ln -sfT /usr/include/hdf5/serial include
ln -sfT /usr/lib/x86_64-linux-gnu/hdf5/serial lib
popd

ls -l /usr/local/hdf5

cd H5Part-1.6.1

./configure --prefix=$PREFIX/h5part --with-hdf5path=/usr/local/hdf5

make -j8 install

cd ..

wget http://download.savannah.gnu.org/releases/named-constant/named-constant-0.3.1.tar.gz

tar -xzf named-constant-0.3.1.tar.gz

cd named-constant-0.3.1
make prefix=$PREFIX install

cd ..

export GENNC_HOME=$PREFIX/gennc
export PATH=$PATH:$GENNC_HOME/bin

export CXXFLAGS="-Wall -Wextra -m64 -O3 -march=native"
export CPPFLAGS="-DNDEBUG"
export LDFLAGS="-m64"

./bootstrap

./configure --prefix=$PREFIX/sph-v1.0.0 --with-h5part=$PREFIX/h5part --with-hdf5=/usr/local/hdf5

make -j8 install

export PATH=$PATH:$PREFIX/sph-v1.0.0/bin

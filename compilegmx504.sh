export CXX=PATH_TO/mpicxx
export CC=PATH_TO/mpicc

cd build/
$HOME/bin/cmake/bin/cmake ../ \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DREGRESSIONTEST_DOWNLOAD=OFF \
    -DREGRESSIONTEST_PATH=../regressiontests-5.0.4 \
    -DCMAKE_INSTALL_PREFIX=$HOME/bin/gmxtmp \
    -DGMX_MPI=ON \
    -DGMX_BUILD_MDRUN_ONLY=OFF \
    -DGMX_DOUBLE=OFF \
    -DGMX_OPENMP=ON \
    -DGMX_DEFAULT_SUFFIX=OFF \ 

make #-j 8
#make install

exit 0


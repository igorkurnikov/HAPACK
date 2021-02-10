# build.sh for test-package
# Sean Greenslade, July 2020

set -e

echo Running build...
mkdir -p ${PREFIX}/bin
mkdir -p ${PREFIX}/lib/python3/site-packages/interx
mkdir -p ${PREFIX}/opt/openmpi
mkdir -p ${PREFIX}/opt/interx/test-scripts
mkdir -p ${PREFIX}/opt/interx/FFDB/760
mkdir -p ${PREFIX}/opt/interx/FFDB/890
mkdir -p ${PREFIX}/opt/interx/FFDB/ALL
mkdir -p ${PREFIX}/opt/interx/doc
mkdir -p ${PREFIX}/opt/interx/MATLAB
mkdir -p ${PREFIX}/opt/interx/HIN

# Build TypifyQMPFF
cd typQMPFF
export SVN_REVISION=$(svn info | grep "Revision" | awk '{print $2}')
./c-make-run.sh Cluster2 FloatDoubleFixed on on TypifyQMPFF off
cp ./Utilities/TypifyQMPFF/TypifyQMPFF ${PREFIX}/bin/
cd ..

# TypifyQMPFF rules file.
cp params_typ/TypifyQMPFF.rul ${PREFIX}/opt/interx/

ARBNAME1=ArbalestCL2-FloatDoubleFixed-Cuda-MPI-march-sandybridge.r3221
ARBNAME2=ArbalestSoftening-FloatDoubleFixed-Cuda-MPI-march-sandybridge.r3302

echo Copying Arbalest binaries...
cp -a ${ARBNAME1} ${PREFIX}/bin/
cp -a ${ARBNAME2} ${PREFIX}/bin/
strip ${PREFIX}/bin/${ARBNAME1}
strip ${PREFIX}/bin/${ARBNAME2}
echo Copying Arbalest-latest symlink...
ln -s ${PREFIX}/bin/${ARBNAME1} ${PREFIX}/bin/Arbalest-latest

echo Copying OpenMPI...
tar -xvzf openmpi.tgz -C ${PREFIX}/opt/

echo Copying InterX python libs...
cp -a pylibs/* ${PREFIX}/lib/python3/site-packages/interx/

echo Copying documentation...
cp -a doc/* ${PREFIX}/opt/interx/doc/

echo Copying matlab scripts...
cp -a MATLAB/* ${PREFIX}/opt/interx/MATLAB/

echo Copying helper scripts...
cp -a helper_py/* ${PREFIX}/bin/
cp -a conda_environment.sh ${PREFIX}/bin/

echo Copying Parameters DBs...
cp FFConfig.xml ${PREFIX}/opt/interx/FFDB/
cp -a FFDB-760/* ${PREFIX}/opt/interx/FFDB/760/
cp -a FFDB-890/* ${PREFIX}/opt/interx/FFDB/890/
cp -a FFDB-ALL/* ${PREFIX}/opt/interx/FFDB/ALL/

echo Copying HINs...
tar -xvzf HINs_2020-12-18.tar.gz -C ${PREFIX}/opt/interx/HIN

echo Copying test scripts...
cp -ar test-scripts/* ${PREFIX}/opt/interx/test-scripts/

exit 0

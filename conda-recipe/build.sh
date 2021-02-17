# build.sh HARLEM linux
#set -e

echo Running build...
mkdir -p ${PREFIX}/bin
mkdir -p ${SP_DIR}/site-packages/molset
mkdir -p ${PREFIX}/opt/harlem
mkdir -p ${PREFIX}/opt/harlem/residues_db
mkdir -p ${PREFIX}/opt/harlem/examples

echo ${SP_DIR}

echo Copying molset library:
#cp -rp ${RECIPE_DIR}/../HARLEM/molset/*  ${SP_DIR}/molset
cp -rp ${RECIPE_DIR}/../../BUILD_HARLEM/HARLEMLL/.libs/lib_molsetc.so  ${SP_DIR}/molset/_molsetc.so
cp -rp ${RECIPE_DIR}/../../BUILD_HARLEM/HARLEMLL/molsetc.py  ${SP_DIR}/molset/

cp -a ${RECIPE_DIR}/../examples/*     ${PREFIX}/opt/harlem/examples
cp -a ${RECIPE_DIR}/../residues_db/*  ${PREFIX}/opt/harlem/residues_db

exit 0

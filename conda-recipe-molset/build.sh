# build.sh HARLEM/MOLSET linux
set -x -e

os_id=`lsb_release -i` 
os_id=${os_id:16}
os_id=`echo "$os_id" | xargs`

echo "SITE_PACKAGE_DIR= ", ${SP_DIR}
echo "OS_ID= ", $os_id
echo "PY_VER= ",$PY_VER 

echo Running build...
mkdir -p ${PREFIX}/bin
mkdir -p ${SP_DIR}/molset
mkdir -p ${SP_DIR}/wx
mkdir -p ${PREFIX}/opt
mkdir -p ${PREFIX}/opt/harlem
mkdir -p ${PREFIX}/opt/harlem/residues_db
mkdir -p ${PREFIX}/opt/harlem/examples
mkdir -p ${PREFIX}/opt/harlem/basis

echo "Copying molset and wx libraries:"
if [[ "$PY_VER" = "3.6" ]] && [[ "$os_id"  = "CentOS" ]]; then
  echo "CENTOS OS :  PY_VER = 3.6"
  cp -rp ${RECIPE_DIR}/../../wxPython-4.1.1_GCC_4.8_PY36/wx/*    ${SP_DIR}/wx
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_CENTOS_7.7/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_CENTOS_7.7/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.7" ]] && [[ "$os_id"  = "CentOS" ]]; then
  echo "CENTOS OS :  PY_VER = 3.7"
  cp -rp ${RECIPE_DIR}/../../wxPython-4.1.1_GCC_4.8_PY37/wx/*    ${SP_DIR}/wx
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_CENTOS_7.7/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_CENTOS_7.7/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.8" ]] && [[ "$os_id"  = "CentOS" ]]; then
  echo "CENTOS OS :  PY_VER = 3.8"
  cp -rp ${RECIPE_DIR}/../../wxPython-4.1.1_GCC_4.8_PY38/wx/*    ${SP_DIR}/wx
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_CENTOS_7.7/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_CENTOS_7.7/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.9" ]] && [[ "$os_id"  = "CentOS" ]]; then
  echo "CENTOS OS :  PY_VER = 3.8"
  cp -rp ${RECIPE_DIR}/../../wxPython-4.1.1_GCC_12.1_PY39/wx/*    ${SP_DIR}/wx
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_CENTOS_7.7/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_CENTOS_7.7/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.9" ]] && [[ "$os_id"  = "Ubuntu" ]]; then
  echo "UBUNTU OS :  PY_VER = 3.9"
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_UBUNTU_20.04/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_UBUNTU_20.04/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.8" ]] && [[ "$os_id"  = "Ubuntu" ]]; then
  echo "UBUNTU OS :  PY_VER = 3.8"
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_UBUNTU_20.04/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_UBUNTU_20.04/*  ${PREFIX}/bin
elif [[ "$PY_VER" = "3.7" ]] && [[ "$os_id"  = "Ubuntu" ]]; then 
  echo "UBUNTU OS :  PY_VER = 3.7"
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_UBUNTU_18.04/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_UBUNTU_18.04/*  ${PREFIX}/bin
else
  echo "OS UNKNOWN : PY_VER = UNKNOWN "
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/LIBS_UBUNTU_20.04/*  ${PREFIX}/lib
  cp -rp ${RECIPE_DIR}/../../MOLSET_BIN_LIBS/BIN_UBUNTU_20.04/*  ${PREFIX}/bin
fi

cp -rp ${RECIPE_DIR}/../HARLEM/molset/*  ${SP_DIR}/molset
cp -L ${RECIPE_DIR}/../../BUILD_HARLEM/HARLEMLL/.libs/lib_molsetc.so  ${SP_DIR}/molset/_molsetc.so
cp -rp ${RECIPE_DIR}/../../BUILD_HARLEM/HARLEMLL/molsetc.py  ${SP_DIR}/molset/
cp -p ${RECIPE_DIR}/../HARLEM/linux/harlem_conda   ${PREFIX}/bin/harlem
cp -a ${RECIPE_DIR}/../examples/*     ${PREFIX}/opt/harlem/examples
cp -a ${RECIPE_DIR}/../basis/*        ${PREFIX}/opt/harlem/basis
cp -a ${RECIPE_DIR}/../residues_db/*  ${PREFIX}/opt/harlem/residues_db

if [[ "$PY_VER" = "3.9" ]] && [[ "$os_id"  = "Ubuntu" ]]; then 
  pip install wxPython-4.1.1-cp39-cp39-linux_x86_64.whl
elif [[ "$PY_VER" = "3.8" ]] && [[ "$os_id"  = "Ubuntu" ]]; then 
  pip install wxPython-4.1.1-cp38-cp38-linux_x86_64.whl
elif [[ "$PY_VER" = "3.7" ]] && [[ "$os_id"  = "Ubuntu" ]]; then 
  pip install wxPython-4.1.1-cp37-cp37m-linux_x86_64.whl
else
  echo "Do not install wxPython wheel" 
fi
exit 0

#!/bin/sh

/bin/rm -f leap.in agl_solv_sleap.prmtop agl_solv_sleap.xyz

cat > leap.in << EOF
set default echo on
source leaprc.amoeba
agl = loadpdb agl.pdb
solvatebox agl WATBOX 6.0
savepdb agl agl_solv.pdb
saveamoebaparm agl agl_solv_sleap.prmtop agl_solv_sleap.crd
quit
EOF

sleap11 < leap.in > agl_solv_sleap.out


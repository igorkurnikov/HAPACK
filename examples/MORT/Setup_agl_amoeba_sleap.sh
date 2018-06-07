#!/bin/sh

/bin/rm -f leap.in agl_sleap.prmtop agl_sleap.xyz

cat > leap.in << EOF
set default echo on
source leaprc.amoeba
agl = loadpdb agl.pdb
saveamoebaparm agl agl_sleap.prmtop agl_sleap.xyz
quit
EOF

sleap11 < leap.in > agl_sleap.out


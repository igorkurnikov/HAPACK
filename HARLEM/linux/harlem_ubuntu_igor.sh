#!/bin/sh
export HARLEM_HOME="/share/apps/HARLEM"
export PYTHONPATH="$HARLEM_HOME:$PYTHONPATH"
python3 -m harlempy $@

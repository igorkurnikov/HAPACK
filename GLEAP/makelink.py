#!/usr/bin/env python

from string import split, strip


f = open( "filelist", 'r' )

newdirs = set()
line = f.readline()
while len(line)>0:

    line  = strip( line )
    items = split( line, '/' )



    assert( len(items) >= 4 )
    assert items[0]==""
    assert items[1]=="Users"
    assert items[2]=="zweig"
    assert items[3]=="Development"

    for i in xrange(4, len(items)-1 ):
        if i==4:
            curtdir = items[i]
        else:
            curtdir = curtdir + "/" + items[i]

        if not(curtdir in newdirs):
            print "mkdir " + curtdir
            newdirs.add( curtdir )

    cmd = "cp " + line + " " + curtdir + "/" + items[-1]
    print cmd



    line = f.readline()

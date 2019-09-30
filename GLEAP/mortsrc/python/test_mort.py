#!/usr/bin/env python

from pymort import *

print "testing hash"
assert( ATOM == hash('atom') )

print "\ntesting unhash"
assert( 'atom' == unhash(ATOM) )

print "\ntesting numvec"
v = numvec(3)
v[0] = 1.0
v[1] = 1.0
v[2] = 1.0
print v
for d in v:
    print d
    assert( d == 1.0 )

print "\ntesting entity"
e = entity_t()
e.name = 'earl'
print 'name:', e.name
assert( e.name == 'earl' )

e.id = 0
print 'id:  ', e.id
assert( e.id == 0 )

e.position = v
print 'pos: ', e.position

print "\ntesting molecule and pointer"
m = molecule_t()
a0 = m.create_atom()
a1 = m.create_atom()
a2 = m.create_atom()
b0 = create_bond(a0, a1)
b1 = create_bond(a0, a2)


print 'a0.id:', a0.id
print 'a1.id:', a1.id
print 'a2.id:', a2.id
print 'b0.id:', b0.id
print 'b1.id:', b1.id

a0.name = 'O1'
a1.name = 'H2'
a2.name = 'H3'

a0.type = 'OW'
a1.type = 'HW'
a2.type = 'HW'

b0.order = 1
b1.order = 1

print 'a0.name:', a0.name
print 'a1.name:', a1.name
print 'a2.name:', a2.name

print 'a0.type:', a0.type
print 'a1.type:', a1.type
print 'a2.type:', a2.type

print 'a0.element:', a0.element
print 'a1.element:', a1.element
print 'a2.element:', a2.element

print 'b0.order:', b0.order
print 'b1.order:', b1.order

print 'a0.is_connected_to(a1):', a0.is_connected_to(a1)
print 'a0.is_connected_to(a2):', a0.is_connected_to(a2)
print 'a1.is_connected_to(a2):', a1.is_connected_to(a2)
print 'b0.is_connected_to(a0):', b0.is_connected_to(a0)
print 'b0.is_connected_to(a1):', b0.is_connected_to(a1)
print 'b0.is_connected_to(a2):', b0.is_connected_to(a2)
print 'b1.is_connected_to(a0):', b1.is_connected_to(a0)
print 'b1.is_connected_to(a1):', b1.is_connected_to(a1)
print 'b1.is_connected_to(a2):', b1.is_connected_to(a2)


r0 = m.create_resd()
r0.name = "HOH"
r0.connect(a0)
r0.connect(a1)
r0.connect(a2)
print 'r0.related_atoms.size():', r0.related_atoms.size()

print "m.natom:", m.natom
print "m.nbond:", m.nbond
print "\ntesting range"

print "m.atoms.size():", m.atoms.size()
print "m.bonds.size():", m.bonds.size()

print "m.atoms[0].name:", m.atoms[0].name
print "m.atoms[1].name:", m.atoms[1].name
print "m.atoms[2].name:", m.atoms[2].name

print "m.bonds[0].order:", m.bonds[0].order
print "m.bonds[1].order:", m.bonds[1].order

print "\ntesting for-in loop with range"

for a in m.atoms:
    print "m.atoms[%d] name: %4s" % (a.id-1, a.name)


import pygtk
pygtk.require('2.0')

import gtk
from pymort import *

class action:
    def __init__(self, curtmol, glarea, level, parm):
        self.glarea   = glarea
	self.curtmol  = curtmol
        self.argument = level + " " + parm
        self.graphnme = "label_" + curtmol.name

    def execute(self):
        if self.glarea.has(self.graphnme) :
	    self.glarea.remove(self.graphnme)
	self.glarea.add("label", self.curtmol, self.argument)


class dialog(gtk.Dialog):
    def __init__(self, parent):
        gtk.Dialog.__init__(self, "Welcome to label...", parent,
	                    gtk.DIALOG_MODAL|gtk.DIALOG_DESTROY_WITH_PARENT,
			    (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
			     gtk.STOCK_OK,     gtk.RESPONSE_OK) )

        self.curtmol = parent.current
	self.glarea  = parent.glarea
        self.level   = "atom"
	self.parm    = "name"


	levelabel = gtk.Label("\nWhich level of molecule?")
	self.vbox.pack_start(levelabel, True, True, 0)
	levelabel.show()

        atombutton = gtk.RadioButton(None, "Atom")
	atombutton.connect( "toggled", self.level_change, "Atom" )
        self.vbox.pack_start(atombutton, True, True, 0)
        atombutton.show()

	resdbutton = gtk.RadioButton(atombutton, "Residue")
	resdbutton.connect( "toggled", self.level_change, "Resd" )
	self.vbox.pack_start(resdbutton, True, True, 0)
        resdbutton.show()

	molebutton = gtk.RadioButton(atombutton, "Molecule")
	molebutton.connect( "toggled", self.level_change, "Mol" )
	self.vbox.pack_start(molebutton, True, True, 0)
        molebutton.show()

        separator1 = gtk.HSeparator()
	self.vbox.pack_start(separator1, True, True, 0)
	separator1.show()

        parmlabel = gtk.Label("\nWhich parameter for label?")
	self.vbox.pack_start(parmlabel, True, True, 0)
	parmlabel.show()

        self.parmopt = gtk.combo_box_new_text()
        self.parmopt.append_text( "name" )
	self.parmopt.append_text( "type" )
	self.parmopt.append_text( "pchg" )
	self.parmopt.connect( "changed", self.parm_change )
	self.vbox.pack_start(self.parmopt, True, True, 0)
        self.parmopt.show()
 
        nulllabel = gtk.Label(" ")
	self.vbox.pack_start(nulllabel, True, True, 0)
	nulllabel.show()

    def level_change(self, widget, level):
        self.level = level

    def parm_change(self, widget):
        model  = self.parmopt.get_model()
	active = self.parmopt.get_active()

	if active < 0:
	    return

	self.parm = model[active][0]

    def get_action(self):
        return action(self.curtmol, self.glarea, self.level, self.parm)

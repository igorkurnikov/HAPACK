import pygtk
pygtk.require( '2.0')

import gtk
from pymort import *

class dialog(gtk.Dialog) :
    def __init__(self, parent, mdb) :
        gtk.Dialog.__init__(self, "Welcome to solvate...", parent,
	                    gtk.DIALOG_MODAL|gtk.DIALOG_DESTROY_WITH_PARENT,
			    (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
 			     gtk.STOCK_OK,     gtk.RESPONSE_OK) )
        self.mdb = mdb
        self.mol = parent.current

        typelabel = gtk.Label("\nWhich type of solvation?")
	self.vbox.pack_start(typelabel, True, True, 0)
	typelabel.show()

        boxbutton = gtk.RadioButton(None, "Box")
	boxbutton.connect( "toggled", self.type_change, "Box" )
        self.vbox.pack_start(boxbutton, True, True, 0)
        boxbutton.show()

	capbutton = gtk.RadioButton(boxbutton, "Cap")
	capbutton.connect( "toggled", self.type_change, "Cap" )
	self.vbox.pack_start(capbutton, True, True, 0)
        capbutton.show()

	octbutton = gtk.RadioButton(boxbutton, "Oct")
	octbutton.connect( "toggled", self.type_change, "Oct" )
	self.vbox.pack_start(octbutton, True, True, 0)
        octbutton.show()

        shellbutton = gtk.RadioButton(boxbutton, "Shell")
	shellbutton.connect( "toggled", self.type_change, "Shell")
	self.vbox.pack_start(shellbutton, True, True, 0)
	shellbutton.show()

	seperator1 = gtk.HSeparator()
	self.vbox.pack_start(seperator1, True, True, 0)
	seperator1.show()

        solventlabel = gtk.Label("\nWhat kind of solvent?")
	self.vbox.pack_start(solventlabel, True, True, 0)
	solventlabel.show()

        self.solventopt = gtk.combo_box_new_text()
        for solvent in mdb.keys():
            self.solventopt.append_text(solvent)
	self.solventopt.connect( "changed", self.solvent_change )
	self.vbox.pack_start(self.solventopt, True, True, 0)
        self.solventopt.show()
       
        nulllabel = gtk.Label(" ")
	self.vbox.pack_start(nulllabel, True, True, 0)
	nulllabel.show()

        separator2 = gtk.HSeparator()
	self.vbox.pack_start(separator2, True, True, 0)
	separator2.show()

        self.centerlabel = gtk.Label( "Center position:" )
	self.vbox.pack_start(self.centerlabel, True, True, 0)
	
	self.centermask  = gtk.Entry()
	self.centermask.set_editable(True)
	self.vbox.pack_start(self.centermask, True, True, 0)

	self.centersepa  = gtk.HSeparator()
	self.vbox.pack_start(self.centersepa)

        self.sizelabel = gtk.Label( "\nBox edge:" );
	self.vbox.pack_start(self.sizelabel, True, True, 0)
	self.sizelabel.show()

	self.sizeentry = gtk.Entry( )
        self.sizeentry.set_editable(True)
	self.vbox.pack_end(self.sizeentry, True, True, 0)
        self.sizeentry.show()

    def type_change(self, widget, type):
        self.solventype = type

	if type == 'Cap' :
	    self.centerlabel.show()
	    self.centermask.show()
	    self.centersepa.show()
	else:
	    self.centerlabel.hide()
	    self.centermask.hide()
	    self.centersepa.hide()

    def solvent_change(self, widget) :
        model  = self.solventopt.get_model()
	active = self.solventopt.get_active()
	if active < 0 :
	    return
      
        self.solvent = self.mdb.get_mol( model[active][0] )
	

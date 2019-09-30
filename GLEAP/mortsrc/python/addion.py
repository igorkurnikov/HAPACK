import pygtk
pygtk.require('2.0')

import gtk
from pymort import *

class action:
    def __init__(self, mol, ion, num) :
        self.mol = mol
	self.ion = ion
	self.num = num
    
    def execute(self):
        shell_extent = 4.0
	resolution   = 1.0
        add_ion(self.mol, self.ion, int(self.num), shell_extent, resolution)
        pass

    def undo(self):
        pass

class dialog(gtk.Dialog) :
    def __init__(self, parent, mdb) :
	gtk.Dialog.__init__(self, "Welcome to add ion...", parent,
	                    gtk.DIALOG_MODAL|gtk.DIALOG_DESTROY_WITH_PARENT,
			    (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
 			     gtk.STOCK_OK,     gtk.RESPONSE_OK) )
        self.mdb = mdb
        self.mol = parent.current

        self.label1 = gtk.Label( "\nNow, total charge is: %5.2f\n" % charge(self.mol) );
	self.vbox.pack_start(self.label1, True, True, 0)
	self.label1.show()
 
        separator1 = gtk.HSeparator()
	self.vbox.pack_start(separator1, True, True, 0)
	separator1.show()


        typelabel = gtk.Label("\nWhat kind of ions?")
	self.vbox.pack_start(typelabel, True, True, 0)
	typelabel.show()

        self.opt = gtk.combo_box_new_text()
        for ioname in mdb.keys():
            self.opt.append_text(ioname)
	self.opt.connect( "changed", self.update )
	self.vbox.pack_start(self.opt, True, True, 0)
        self.opt.show()

        numblabel = gtk.Label("\nHow many of them?")
	self.vbox.pack_start(numblabel, True, True, 0)
	numblabel.show()

        adj = gtk.Adjustment(0.0, 0.0, 1000.0, 1.0, 5.0, 0.0)
        self.spinner = gtk.SpinButton(adj, 0, 0)
	self.spinner.connect( "changed", self.update )
        self.vbox.pack_start(self.spinner, False, True, 0)
	self.spinner.show()

        space = gtk.Label("\n")
	self.vbox.pack_start(space,True,True,0)
        space.show()

	separator2 = gtk.HSeparator()
	self.vbox.pack_start(separator2, True, True, 0)
	separator2.show()

        self.label2 = gtk.Label( "\nTotal charge will be: %5.2f\n" % charge(self.mol)  )
        self.vbox.pack_end(self.label2, True, True, 0)
        self.label2.show()

	self.vbox.show()

    def update(self, wdiget):
        model = self.opt.get_model()
	active = self.opt.get_active()
	if active < 0 :
	    return
      
        self.ion = self.mdb.get_mol( model[active][0] )
	self.num = self.spinner.get_value()
        totpchg = charge(self.mol) + self.num * charge(self.ion)
	self.label2.set_text( "\nTotal charge will be: %5.2f\n" % totpchg )


    def get_action(self):
        return action( self.mol, self.ion, self.num )



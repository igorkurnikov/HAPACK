import pygtk
pygtk.require('2.0')

import gtk
from pymort import *


class action:
    def __init__(self, curtmol, glarea, graftype, oper, part):
        self.grafname = graftype + "_" + curtmol.name
	self.curtmol = curtmol
	self.glarea = glarea
	self.oper = oper
	self.part = part

    def execute(self):
        graf = self.glarea.get( self.grafname )
	if graf is None:
            return

	assert self.oper=="on" or self.oper=="off" or self.oper=="only"

	if self.oper == "only":
            graf.setall(VISIBLE, OFF)
    
	atms = mask_atom(self.curtmol, self.part)

	if self.oper == "only" or self.oper == "on":
	    graf.set(atms, VISIBLE, ON)
	else:
	    graf.set(atms, VISIBLE, OFF)


class dialog(gtk.Dialog):
    def __init__(self, parent):
        gtk.Dialog.__init__(self, "Welcome to show/hide...", parent,
	                    gtk.DIALOG_MODAL|gtk.DIALOG_DESTROY_WITH_PARENT,
			    (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
			     gtk.STOCK_OK,     gtk.RESPONSE_OK) )

        self.curtmol = parent.current
	self.glarea = parent.glarea
        self.graf = "molecule"
	self.oper = "on"
	self.part = "alias.all"
	
        graflabel = gtk.Label("\nWhich graphic to show/hide?")
	self.vbox.pack_start(graflabel, True, True, 0)
	graflabel.show()

        selfbutton = gtk.RadioButton(None, "Itself")
	selfbutton.connect( "toggled", self.graf_change, "molecule" )
	self.vbox.pack_start(selfbutton, True, True, 0)
	selfbutton.show()

        labelbutton = gtk.RadioButton(selfbutton, "Label")
	labelbutton.connect( "toggled", self.graf_change, "label" )
        self.vbox.pack_start(labelbutton, True, True, 0)
	labelbutton.show()

        separator1 = gtk.HSeparator()
	self.vbox.pack_start(separator1, True, True, 0)
	separator1.show()

        operlabel = gtk.Label("\nWhat operation to perform?")
	self.vbox.pack_start(operlabel, True, True, 0)
	operlabel.show()

        onbutton = gtk.RadioButton(None, "on")
	onbutton.connect( "toggled", self.oper_change, "on" )
	self.vbox.pack_start(onbutton, True, True, 0)
	onbutton.show()

        offbutton = gtk.RadioButton(onbutton, "off")
	offbutton.connect( "toggled", self.oper_change, "off" )
	self.vbox.pack_start(offbutton, True, True, 0)
	offbutton.show()

        onlybutton = gtk.RadioButton(onbutton, "only")
	onlybutton.connect( "toggled", self.oper_change, "only" )
	self.vbox.pack_start(onlybutton, True, True, 0)
	onlybutton.show()

        separator2 = gtk.HSeparator()
	self.vbox.pack_start(separator2, True, True, 0)
	separator2.show()

        partlabel = gtk.Label("\nWhich part to show/hide?")
        self.vbox.pack_start(partlabel, True, True, 0)
	partlabel.show()

	self.partopt = gtk.combo_box_new_text()
        self.partopt.append_text( "all" )
	self.partopt.append_text( "hydrogens" )
	self.partopt.append_text( "backbone" )
	self.partopt.append_text( "Ambmask..." )
	self.partopt.append_text( "Smiles...." )
        self.partopt.connect( "changed", self.part_change )
	self.vbox.pack_start(self.partopt, True, True, 0)
        self.partopt.show()

	self.masktxt = gtk.Entry()
        self.vbox.pack_start(self.masktxt, True, True, 0)
	
	nulllabel = gtk.Label(" ")
	self.vbox.pack_start(nulllabel, True, True, 0)
	nulllabel.show()

    def graf_change(self, widget, graf):
        self.graf = graf

    def oper_change(self, widget, oper):
        self.oper = oper

    def part_change(self, widget):
        model = self.partopt.get_model()
	active = self.partopt.get_active()

	if active<0:
	    return

	if model[active][0] == "Ambmask..." or model[active][0] == "Smiles....":
	    self.masktxt.show()
	else:
	    self.masktxt.hide()
    
    def get_action(self):
        model = self.partopt.get_model()
	active = self.partopt.get_active()

	if active<0:
	    return None

	if model[active][0] == "Ambmask...": 
	    part = "ambmask." + self.masktxt.get_text()
	elif model[active][0] == "Smiles....":
            part = "smiles." + self.masktxt.get_text()
	else:
	    part = "alias." + model[active][0] 

        return action(self.curtmol, self.glarea, self.graf, self.oper, part)

        

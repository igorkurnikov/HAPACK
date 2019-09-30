#!/usr/bin/env python

import pygtk
pygtk.require('2.0')

import gtk
import gtk.gtkgl

from pymort import *
from OpenGL.GL import *

ui_string = """<ui>
  <menubar name='Menubar'>
    <menu action='FileMenu'>
      <menuitem action='New'/>
      <menuitem action='Open'/>
      <menuitem action='Save'/>
      <separator/>
      <menuitem action='Close'/>
      <menuitem action='Quit'/>
    </menu>
    <menu action='EditMenu'>
      <menuitem action='Addion'/>
      <menuitem action='Solvate'/>
    </menu>
    <menu action='ViewMenu'>
      <menuitem action='Show'/>
      <menuitem action='Label'/>
      <menuitem action='Color'/>
    </menu>
    <menu action='HelpMenu'>
      <menuitem action='About'/>
    </menu>
  </menubar>
  <toolbar name='Toolbar'>
    <toolitem action='New'/>
    <toolitem action='Open'/>
    <toolitem action='Save'/>
    <separator/>
    <toolitem action='Quit'/>
  </toolbar>
</ui>"""

def get_molname(fn):
    try:
        start = fn.rindex('/')+1
    except ValueError:
        start = 0

    end   = fn.rindex('.')
    return fn[start:fn.index('.')]


class glarea_t(gtk.gtkgl.DrawingArea):
    def __init__(self):
        try:
            # try double-buffered
            glconfig = gtk.gdkgl.Config(mode=(gtk.gdkgl.MODE_RGB    |
	                                      gtk.gdkgl.MODE_DOUBLE |
					      gtk.gdkgl.MODE_DEPTH))
        except gtk.gdkgl.NoMatches:
            # try single-buffered
            glconfig = gtk.gdkgl.Config(mode=(gtk.gdkgl.MODE_RGB    |
	                                      gtk.gdkgl.MODE_DEPTH))

	print "glconfig.is_rgba() =",            glconfig.is_rgba()
	print "glconfig.is_double_buffered() =", glconfig.is_double_buffered()
	print "glconfig.has_depth_buffer() =",   glconfig.has_depth_buffer()
  	# get_attrib()
	print "gtk.gdkgl.RGBA = %d"         % glconfig.get_attrib(gtk.gdkgl.RGBA)
	print "gtk.gdkgl.DOUBLEBUFFER = %d" % glconfig.get_attrib(gtk.gdkgl.DOUBLEBUFFER)
	print "gtk.gdkgl.DEPTH_SIZE = %d"   % glconfig.get_attrib(gtk.gdkgl.DEPTH_SIZE)

        gtk.gtkgl.DrawingArea.__init__(self,glconfig)

        #self.add_events(gtk.gdk.BUTTON_PRESS   |
	#                gtk.gdk.MOTION_NOTIFY  |
        #                gtk.gdk.BUTTON_RELEASE )
        self.add_events(gtk.gdk.BUTTON_PRESS_MASK |
	                gtk.gdk.POINTER_MOTION_MASK |
			gtk.gdk.BUTTON_RELEASE_MASK)

        self.connect_after('realize', self.initgl)
        self.connect('configure_event', self.resize)
        self.connect('expose_event', self.paint)
	self.connect('button_press_event', self.mouse_press)
	self.connect('button_release_event', self.mouse_release)
	self.connect('motion_notify_event', self.mouse_move)

        #self.connect('map_event', self.map)
        self.impl = gldrawing_t()

    def glbegin(self):
        # get GLContext and GLDrawable
	glcontext = self.get_gl_context()
	gldrawable = self.get_gl_drawable()

	# GL calls
	if not gldrawable.gl_begin(glcontext): 
	    print "Can't begin GL"
	    return

    def glend(self):
	gldrawable = self.get_gl_drawable()

	if gldrawable.is_double_buffered():
	    gldrawable.swap_buffers()
	else:
	    glFlush()

        gldrawable.gl_end()


    def initgl(self,widget):
        self.glbegin()
	self.impl.init()
	self.glend()

    def resize(self,widget,event):
        self.glbegin()
	x,y,width,height = self.get_allocation()
        self.impl.resize(width, height)
	self.glend()

    def paint(self,widget,event):
        self.glbegin()
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
	self.impl.repaint()
	self.glend()

    def add(self,type,mol, argument):
        self.impl.add(type, mol, argument)
        self.paint(self,None)

    def has(self, name):
        return self.impl.has(name)

    def get(self, name):
        return self.impl.get(name)

    def remove(self, name):
        self.impl.remove(name)
    
    def mouse_press(self,widget,event):
        self.impl.mouse_press(event.button, event.state, event.x, event.y)

    def mouse_move(self,widget,event):
        #self.impl.mouse_move(0, event.state, event.x, event.y)
        #self.paint(self,None)
        pass

    def mouse_release(self,widget,event):
        self.impl.mouse_release(event.button, event.state, event.x, event.y)
        self.paint(self,None)

class mainwin_t(gtk.Window):

    def __init__(self):
        gtk.Window.__init__(self)
	self.set_position(gtk.WIN_POS_CENTER)
	self.set_title('PyLeap')
        self.connect('delete_event', self.delete)
        self.set_size_request(640,480)
	vbox = gtk.VBox()
	self.add(vbox)

	self.create_ui()
	vbox.pack_start(self.ui.get_widget('/Menubar'), expand=False)
	vbox.pack_start(self.ui.get_widget('/Toolbar'), expand=False)

        self.glarea = glarea_t()
	vbox.pack_start(self.glarea)

	status = gtk.Statusbar()
	vbox.pack_end(status, expand=False)

        self.current = None
        self.history = []

    def initdb(self, amberhome)	:
        self.content = database_t()
        pgwin = gtk.Window(gtk.WINDOW_TOPLEVEL)
	pgwin.set_size_request(300,100)
        vbox = gtk.VBox(False, 5)
   	vbox.set_border_width(10)
     	pgwin.add(vbox)
        vbox.show()
	    	  
	# Create label
	info1 = gtk.Label("Initilzating:")
	vbox.pack_start(info1, False, False, 0)
	info1.show()

	info2 = gtk.Label("             ")
	vbox.pack_start(info2, False, False, 0)
        info2.show()

	# Create the ProgressBar
	pgbar = gtk.ProgressBar()
        vbox.pack_end(pgbar)
	pgbar.show()

	pgwin.show_all()


        info2.set_text("loading amino acids templates")
        aminoacids  = load_mdb( amberhome + "/dat/leap/lib/all_amino94.lib")
        self.content.set( "_aminoacids", aminoacids )
	pgbar.set_fraction(0.25)

        info2.set_text("loading C-terminal amino acids templates")
	aminoacidsct = load_mdb(amberhome + "/dat/leap/lib/all_aminoct94.lib")
        self.content.set( "_aminoacidsct", aminoacidsct)
	pgbar.set_fraction(0.50)

        info2.set_text("loading N-terminal amino acids templates")
        aminoacidsnt = load_mdb(amberhome + "/dat/leap/lib/all_aminont94.lib")
	self.content.set( "_aminoacidsnt", aminoacidsnt)
	pgbar.set_fraction(0.75)


        info2.set_text("loading solvents and ions templates")
        solvents = load_mdb( amberhome + "/dat/leap/lib/solvents.lib" )
        ions = load_mdb( amberhome + "/dat/leap/lib/ions94.lib")
	self.content.set( "_solvents", solvents )
	self.content.set( "_ions", ions )
        pgbar.set_fraction(1.00)

        pgwin.destroy()

    def create_ui(self):
        ag = gtk.ActionGroup('WindowActions')
        actions = [
            ('FileMenu', None, '_File'),
            ('New',      gtk.STOCK_NEW, '_New', '<control>N',
             'Create a new file', self.file_new_cb),
            ('Open',     gtk.STOCK_OPEN, '_Open', '<control>O',
             'Open a file', self.file_open_cb),
	    ('Save',     gtk.STOCK_SAVE, '_Save', '<control>S',
	     'Save a file', self.file_save_cb),
            ('Close',    gtk.STOCK_CLOSE, '_Close', '<control>W',
             'Close the current window', self.file_close_cb),
            ('Quit',     gtk.STOCK_QUIT, '_Quit', '<control>Q',
             'Quit application', self.delete),
	    ('EditMenu', None, '_Edit'),
	    ('Addion', None, '_Add ions', None,
	     'Add counter ion', self.edit_addion_cb),
	    ('Solvate', None, '_Solvate', None,
	     'Add solvation',  self.edit_solvate_cb),
	    ('ViewMenu', None, '_View'),
	    ('Show', None, '_Show/Hide', None,
	     'Show', self.view_show_cb),
            ('Label', None, '_Label', None,
	     'Label', self.view_label_cb),
	    ('Color', None, '_Color', None,
	     'Color', self.view_color_cb),
            ('HelpMenu', None, '_Help'),
            ('About',    None, '_About', None, 'About application',
             self.help_about_cb),
            ]
        ag.add_actions(actions)
        self.ui = gtk.UIManager()
        self.ui.insert_action_group(ag, 0)
        self.ui.add_ui_from_string(ui_string)
        self.add_accel_group(self.ui.get_accel_group())

    def file_new_cb(self, action):
        pass

    def file_save_cb(self, action):
        if self.current is None:
	    self.popup_message( "Error: no molecule to save" )
	    return

        dialog = gtk.FileChooserDialog("Save..", self,
	                               gtk.FILE_CHOOSER_ACTION_SAVE,
				       (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
				        gtk.STOCK_SAVE,   gtk.RESPONSE_OK))
        filter = gtk.FileFilter()
        filter.set_name("Molecular files")
        filter.add_pattern("*.sdf")
	filter.add_pattern("*.pdb")
	filter.add_pattern("*.mol2")
        dialog.add_filter(filter)
        dialog.hide()

        if dialog.run() == gtk.RESPONSE_OK :
           fn = dialog.get_filename()
	   save_mol(fn, self.current)

        dialog.destroy()

    def file_open_cb(self, action):
        dialog = gtk.FileChooserDialog("Open..", self,
                                       gtk.FILE_CHOOSER_ACTION_OPEN,
                                       (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                        gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)

        filter = gtk.FileFilter()
        filter.set_name("Molecular files")
        filter.add_pattern("*.sdf")
	filter.add_pattern("*.pdb")
	filter.add_pattern("*.mol2")
        dialog.add_filter(filter)
        dialog.hide()

        if dialog.run() == gtk.RESPONSE_OK:
	    fn = dialog.get_filename()
	    m = load_mol( fn, self.content )
            m.name = get_molname(fn)
	    self.glarea.add("molecule", m, "line")
	    self.content.set(m.name, m)
            self.current = m

        dialog.destroy()

    def file_close_cb(self, action):
        pass

    def edit_addion_cb(self, action):
        import addion
        if self.current is None:
	    self.popup_message("Error: no molecule to add ions to.")
	    return

        dialog = addion.dialog(self, self.content.get_mdb("_ions") )

        if dialog.run()==gtk.RESPONSE_OK:
	   action = dialog.get_action()

           if action.execute():
	       self.history.append( action )
               self.glarea.paint(self.glarea, None)

        dialog.destroy()

    def edit_solvate_cb(self,action):
        import solvate
	if self.current is None:
	    self.popup_message("Error: no molecule to add solvation to.")
	    return

	dialog = solvate.dialog(self, self.content.get_mdb("_solvents") )

	if dialog.run()==gtk.RESPONSE_OK:
	    pass

	dialog.destroy()

    def view_show_cb(self, action):
        if self.current is None:
	    self.popup_message("Error: no moleucle has been loaded")
            return

        import show
	dialog = show.dialog(self)
	if dialog.run()==gtk.RESPONSE_OK:
	    action = dialog.get_action()
	    action.execute()

	dialog.destroy()

    def view_label_cb(self, action):
        if self.current is None:
	    self.popup_message("Error: no moleucle has been loaded")
            return

	import label
	dialog = label.dialog(self)

        if dialog.run()==gtk.RESPONSE_OK:
	   action = dialog.get_action()
	   action.execute()

	dialog.destroy()
        pass

    def view_color_cb(self, action):
        pass

    def popup_message(self, message):
        dialog = gtk.MessageDialog(self,
                                   (gtk.DIALOG_MODAL |
                                    gtk.DIALOG_DESTROY_WITH_PARENT),
                                   gtk.MESSAGE_INFO, gtk.BUTTONS_OK,
                                   message)
        dialog.run()
        dialog.destroy()

    def help_about_cb(self, action):
        self.popup_message( "PyLeap" )

    def delete(self, widget, event=None, data=None):
        dialog = gtk.Dialog("Confirmation", self, 
	                gtk.DIALOG_MODAL|gtk.DIALOG_DESTROY_WITH_PARENT, 
	                (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, 
			 gtk.STOCK_OK,     gtk.RESPONSE_OK) )

        label = gtk.Label("\nReally Quit?\n")
	label.show()
        dialog.vbox.pack_start(label, True, True, 0)
  
	if( dialog.run() == gtk.RESPONSE_OK ) :
	    gtk.main_quit()

        dialog.destroy()



if __name__ == "__main__":
    import os
    if os.getenv('AMBERHOME') is None:
        print "Error: environment AMBERHOME has not been set!"
        from sys import exit
        exit(-1)

    amberhome = os.getenv('AMBERHOME')
    win = mainwin_t()
    win.initdb(amberhome)
    win.show_all()
    gtk.main()


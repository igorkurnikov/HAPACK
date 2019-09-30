static const int DEFAULT_BOX_PADDING_SIZE = 0;


#define ACTION_HINT_FILE_OPEN "Open file"
#define ACTION_HINT_FILE_CLOSE "Close data" 
#define ACTION_HINT_FILE_SAVE "Save file"
#define ACTION_HINT_FILE_QUIT "Quit the application"

#define ACTION_HINT_EDIT_FIXBOND "Fix bond order automatically"
#define ACTION_HINT_EDIT_ADDHYDR "Add hydrogens"
#define ACTION_HINT_EDIT_SETPCHG "Set partial charge"
#define ACTION_HINT_EDIT_PARMCHK "Run parmchk to get force field parameter"

#define ACTION_HINT_EDIT_FASTBLD "Fast build polymer using library and fixbo"
#define ACTION_HINT_EDIT_ADDIONS "Add counter ions"
#define ACTION_HINT_EDIT_SOLVATE "Add solvation"
#define ACTION_HINT_EDIT_FFPARAM "Save force field parameters"

#define ACTION_HINT_VIEW_SHIDE "Show/Hide atoms"
#define ACTION_HINT_VIEW_LABEL "Show/Hide labels"
#define ACTION_HINT_VIEW_COLOR "Set color"
#define ACTION_HINT_VIEW_RIBBON "Show ribbon"
#define ACTION_HINT_VIEW_SURFACE "Show surface"

DEFINE_MENU_ITEM(file_open, ACTION_HINT_FILE_OPEN)
DEFINE_MENU_ITEM(file_close,ACTION_HINT_FILE_CLOSE)
DEFINE_MENU_ITEM(file_save, ACTION_HINT_FILE_SAVE)
DEFINE_MENU_ITEM(file_quit, ACTION_HINT_FILE_QUIT)

DEFINE_MENU_ITEM(edit_fixbond, ACTION_HINT_EDIT_FIXBOND)
DEFINE_MENU_ITEM(edit_addhydr, ACTION_HINT_EDIT_ADDHYDR)
DEFINE_MENU_ITEM(edit_setpchg, ACTION_HINT_EDIT_SETPCHG)
DEFINE_MENU_ITEM(edit_parmchk, ACTION_HINT_EDIT_PARMCHK)

DEFINE_MENU_ITEM(edit_fastbld, ACTION_HINT_EDIT_FFPARAM)
DEFINE_MENU_ITEM(edit_addions, ACTION_HINT_EDIT_ADDIONS)
DEFINE_MENU_ITEM(edit_solvate, ACTION_HINT_EDIT_SOLVATE)
DEFINE_MENU_ITEM(edit_ffparam, ACTION_HINT_EDIT_FFPARAM)


DEFINE_MENU_ITEM(view_shide, ACTION_HINT_VIEW_SHIDE)
DEFINE_MENU_ITEM(view_label, ACTION_HINT_VIEW_LABEL)
DEFINE_MENU_ITEM(view_color, ACTION_HINT_VIEW_COLOR)
DEFINE_MENU_ITEM(view_ribbon, ACTION_HINT_VIEW_RIBBON)
DEFINE_MENU_ITEM(view_surface, ACTION_HINT_VIEW_SURFACE)

GtkActionEntry MAINWIN_ACTION_ENTRYS[] =
{
    { "FileMenuAction", NULL, "_File"},   // name, stockid, label
    { "OpenAction", GTK_STOCK_OPEN, "_Open",//name, stockid, label
      "<control>O", ACTION_HINT_FILE_OPEN,    // accelerate, tooltip
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(file_open))//callback
    },
    { "CloseAction", NULL, "Close",
      NULL, ACTION_HINT_FILE_CLOSE,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(file_close))
    }, 
    { "SaveAction", GTK_STOCK_SAVE, "_Save",
      "<control>S", ACTION_HINT_FILE_SAVE,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(file_save))
    },
    { "QuitAction", GTK_STOCK_QUIT, "_Quit",
      "<control>Q", ACTION_HINT_FILE_QUIT,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(file_quit))
    },
    { "EditMenuAction", NULL, "_Edit"},
    { "MonomerAction", NULL, "Monomer"},
    { "FixbondAction", NULL, "Fix bond",
      NULL, ACTION_HINT_EDIT_FIXBOND,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_fixbond))
    },
    { "AddhydrAction", NULL, "Add hydrogen",
      NULL, ACTION_HINT_EDIT_ADDHYDR,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_addhydr))
    },
    { "SetpchgAction", NULL, "Set charge",
      NULL, ACTION_HINT_EDIT_SETPCHG,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_setpchg))
    },
    { "ParmchkAction", NULL, "Run parmchk",
      NULL, ACTION_HINT_EDIT_PARMCHK,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_parmchk))
    },
    { "PolymerAction", NULL, "Polymer"},
    { "FastbldAction", NULL, "Fast build",
      NULL, ACTION_HINT_EDIT_FASTBLD,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_fastbld))
    },
    { "AddionsAction", NULL, "Add ions",
      NULL, ACTION_HINT_EDIT_ADDIONS,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_addions))
    },
    { "SolvateAction", NULL, "Solvate",
      NULL, ACTION_HINT_EDIT_SOLVATE,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_solvate))
    },
    { "FfparamAction", NULL, "Save ffparm",
      NULL, ACTION_HINT_EDIT_FFPARAM,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(edit_ffparam))
    },
    
    { "ViewMenuAction", NULL, "_View"},
    { "ShideAction", NULL, "Show/Hide",
      NULL, ACTION_HINT_VIEW_SHIDE,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(view_shide))
    },
    { "LabelAction", NULL, "Label",
      NULL, ACTION_HINT_VIEW_LABEL,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(view_label))
    },
    { "ColorAction", NULL, "Color",
      NULL, ACTION_HINT_VIEW_COLOR,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(view_color))
    },
    { "RibbonAction", NULL, "Ribbon",
      NULL, ACTION_HINT_VIEW_RIBBON,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(view_ribbon))
    },
    { "SurfaceAction", NULL, "Surface",
      NULL, ACTION_HINT_VIEW_SURFACE,
      G_CALLBACK(STATIC_ACTION_LISTENER_NAME(view_surface))
    }
    
    
};


static const int MAINWIN_ACTION_NENTRY = G_N_ELEMENTS( MAINWIN_ACTION_ENTRYS );

const gchar* MAINWIN_UIDESIGN =
    "<ui>"
    "   <menubar name='Menubar'>"
    "     <menu action='FileMenuAction'>"
    "       <menuitem action='OpenAction'/>"
    "       <menuitem action='CloseAction'/>"
    "       <menuitem action='SaveAction'/>"
    "       <separator/>"
    "       <menuitem action='QuitAction'/>"
    "     </menu>"
    "     <menu action='EditMenuAction'>"
    "       <menu action='MonomerAction'>"
    "         <menuitem action='FixbondAction'/>"
    "         <menuitem action='AddhydrAction'/>"
    "         <menuitem action='SetpchgAction'/>"
    "         <menuitem action='ParmchkAction'/>"
    "       </menu>"
    "       <menu action='PolymerAction'>"
    "         <menuitem action='FastbldAction'/>"
    "       </menu>"
    "       <menuitem action='AddionsAction'/>"
    "       <menuitem action='SolvateAction'/>"
    "       <menuitem action='FfparamAction'/>"
    "     </menu>"
    "     <menu action='ViewMenuAction'>"
    "       <menuitem action='ShideAction'/>"
    "       <menuitem action='LabelAction'/>"
    "       <menuitem action='RibbonAction'/>"
    "       <menuitem action='SurfaceAction'/>"
    "     </menu>"
    "   </menubar>"
    "   <toolbar name='Toolbar'>"
    "     <toolitem action='OpenAction'/>"
    "     <toolitem action='SaveAction'/>"
    "     <separator/>"
    "     <toolitem action='QuitAction'/>"
    "   </toolbar>"
    "</ui>";


#ifndef GLEAP_LEAPSRC_BUILDUI_HPP
#define GLEAP_LEAPSRC_BUILDUI_HPP

#include <gtk/gtkstock.h>
#include <gtk/gtkitem.h>
#include <mortgl.hpp>

typedef struct {
  const gchar     *label;
  const gchar     *icon_file;
  GCallback  callback;
} GtkToolButtonEntry;

using mort::drawing;

void S_on_select_button_clicked(GtkButton* button, gpointer userdata)
{
    drawing()->mouse_settype( "select" );
}

void S_on_move_button_clicked(GtkButton* button, gpointer userdata)
{
    drawing()->mouse_settype( "translate" );
}

void S_on_rotate_button_clicked(GtkButton* button, gpointer userdata)
{
    drawing()->mouse_settype( "rotate" );
}

void S_on_scale_button_clicked(GtkButton* button, gpointer userdata)
{
    drawing()->mouse_settype( "scale" );
}


void S_on_bond_button_clicked(GtkButton* button, gpointer userdata)
{
    drawing()->mouse_settype( "drawbond" );
}

void S_on_ring_button_clicked(GtkButton* button, gpointer userdata)
{
    // drawing()->set_mouser( "drawring" );
}


GtkToolButtonEntry
edit_tool_entries[NUM_TOOL] =
{
    {
        "Select rectangle region",
	"tool-select.png",
	G_CALLBACK(S_on_select_button_clicked)
    },
    {
        "Move molecule",
	"tool-move.png",
        G_CALLBACK(S_on_move_button_clicked)
    },
    {
        "Rotate molecule",
	"tool-rotate.png",
        G_CALLBACK(S_on_rotate_button_clicked)
    },
    {
        "Rotate molecule",
	"tool-scale.png",
        G_CALLBACK(S_on_scale_button_clicked)
    },
    {
        "Draw bond",
	"tool-bond.png",
	G_CALLBACK(S_on_bond_button_clicked)
    },
    {
        "Draw ring",
	"tool-r6.png",
	G_CALLBACK(S_on_ring_button_clicked)
    }
};


#endif

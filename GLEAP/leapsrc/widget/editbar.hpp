#ifndef GLEAP_LEAPSRC_BUILDER_HPP
#define GLEAP_LEAPSRC_BUILDER_HPP

#include <gtk/gtk.h>
#include <vector>

enum
{
    SELECT_TOOL,
    MOVE_TOOL,
    ROTATE_TOOL,
    SCALE_TOOL,
    BOND_TOOL,
    RING_TOOL,
    NUM_TOOL
};

class editbar_t
{
public:

    editbar_t();

    virtual ~editbar_t();

    void show();

    operator GtkWidget*()
    {
        return m_self;
    }

private:

    GtkWidget* m_self;

    std::vector< GtkWidget* > m_buttons;
};



#endif



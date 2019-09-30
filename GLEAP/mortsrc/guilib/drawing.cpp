#include <boost/shared_ptr.hpp>
#include "graphic.hpp"
#include "drawing.hpp"


namespace mort
{

    drawing_i::drawing_i()
    {
    }               

    drawing_i::~drawing_i( )
    {
    }

    void drawing_i::addone( const shared_ptr< graphic_i >& graphic )
    {
        if( !has( graphic->name() ) )
        {
            m_graphics.push_back( graphic );
        }
    }

    bool drawing_i::has( const string& gname ) const
    {
        vector< shared_ptr< graphic_i > >::const_iterator i = m_graphics.begin();

        for( ; i != m_graphics.end(); ++i )
        {
            if( (*i)->name() == gname )
            {
                return true;
            }
        }

        return false;
    }

    void drawing_i::remove( const string& gname )
    {
        vector< shared_ptr< graphic_i > >::iterator i = m_graphics.begin();

        for( ; i != m_graphics.end(); ++i )
        {
            if( (*i)->name() == gname )
            {
                m_graphics.erase( i );
                return;
            }
        }
    }

    
} // namespace amber
 

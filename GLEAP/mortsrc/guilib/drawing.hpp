#ifndef MOLVIEW_WINDOW_H
#define MOLVIEW_WINDOW_H

#include <boost/shared_ptr.hpp>
#include <common.hpp>

namespace mort
{

    class graphic_i;
    class molecule_t;

    /// \brief the main view class of gleap
    /// \ingroup leaplib
    class drawing_i
    {
    public:

        /// constructor
        drawing_i();

        /// deconstructor
        virtual ~drawing_i();
    

        /// test if graphic exist
        /// \param name name of the graphic
        bool has(const string& name) const;
    
        /// remove certain graphic
        /// \param name name of the graphic
        void remove(const string& name);

        /// add graphic for display 
        virtual void add(const string& type, const molecule_t& mol)=0;

        void addone(const shared_ptr<graphic_i>& g);

        /// handler for mouse pressed
        virtual void mouse_press( int button, int state, double x, double y )=0;

        /// handler for mouse move
        virtual void mouse_move ( int button, int state, double x, double y )=0;

        /// handler for mouse released
        virtual void mouse_release( int button, int state, double x, double y )=0;

        /// OpenGL initializaiton
        virtual void init() = 0;
 
        /// paint all graphics
        virtual void repaint() = 0;
   
        /// adjust size of the display
        virtual void resize( int width, int height ) = 0;

        /// OpenGL preparation
        /// varies on differnent system
        virtual void glbegin() = 0;
        
        /// OpenGL ending, flush or swap buffer
        /// varies on differnent system
        virtual void glend() = 0;

        /// send an expose signal for repaint.
        /// varies on differnent system
        virtual void expose() = 0;

    protected:
        vector< shared_ptr< graphic_i > > m_graphics;
    };

    class null_drawing : public drawing_i
    {
    public:

        virtual void add(const string&, const molecule_t& ) {}

        /// handler for mouse pressed
        virtual void mouse_press(int, int, double, double) {}

        /// handler for mouse move
        virtual void mouse_move (int, int, double, double) {}

        /// handler for mouse released
        virtual void mouse_release(int, int, double, double) {}

        virtual void init() {}

        virtual void repaint() {}

	virtual void resize(int, int) {}

        virtual void glbegin() {}
        
        virtual void glend() {}
        
        virtual void expose() {};
    };  

    typedef shared_ptr< drawing_i > drawing_ptr;

} // namespace mort
 
#endif


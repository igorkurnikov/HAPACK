#ifndef GLEAP_MORTSRC_MORTGL_DRAWING_HPP
#define GLEAP_MORTSRC_MORTGL_DRAWING_HPP


#include <common.hpp>


namespace mort
{
    class mouser_i;

    struct ortho_t;

    class graphic_t;

    class molecule_t;

    class drawing_t
    {
    public:

        /// constructor
        drawing_t();

        /// deconstructor
        virtual ~drawing_t();

        void add(const string& type, molecule_t& mol, const string& argument);

        bool has(const string& name) const;

	shared_ptr<graphic_t> get(const string& name) const;

	void remove(const string& name);

        void update(const string& name);

        void mouse_press( int button, int state, double x, double y );

        void mouse_move ( int button, int state, double x, double y );

        void mouse_release( int button, int state, double x, double y );

        void mouse_settype( const string& mouser );

        void init();

        void update();
 
        void repaint();
    
        void resize( int width, int height );

        ortho_t& ortho() 
        {
            return *mp_ortho;
        }

    private:

        void addone(const shared_ptr<graphic_t>& graphic);

    private:

        vector< shared_ptr<graphic_t> > m_graphics;
        shared_ptr< mouser_i > mp_mouser;
        shared_ptr< ortho_t > mp_ortho;
        quantity_e m_quantity;
        string m_mouser;
    };

    typedef shared_ptr< drawing_t > drawing_ptr;

    drawing_ptr& drawing();

} // namespace mort


#endif

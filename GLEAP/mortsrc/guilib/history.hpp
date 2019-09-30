#ifndef GTKLEAP_HISTORY_HPP
#define GTKLEAP_HISTORY_HPP

#include <common.hpp>

namespace mort
{
    
    template< typename T >
    class history_T
    {
    public:

        history_T() 
        { 
            m_curt = 0; 
        }

        virtual ~history_T() 
        { 
        }
    
        void push( const T& cmd ) 
        {
            m_cmds.push_back( cmd );    
            m_curt = m_cmds.size();
        }
        
        void inc()
        {
            if( m_curt < (int)m_cmds.size() )
            {
                m_curt++;
            }
        }
        
        void dec()
        {
            if( m_curt > 0 )
            {
                m_curt--;
            }
        }
        
        T val() const
        {
            return m_curt < (int)m_cmds.size() ? m_cmds[ m_curt ] : T();
        }
        
    private:

        vector< T >  m_cmds;
        
        int m_curt;
    };
    
} // namespace mort
    
#endif

#ifndef DDUTILS_MORT_ADJACENCY_H
#define DDUTILS_MORT_ADJACENCY_H

#include <vector>

namespace mort
{
    using std::vector;

    /// \brief ajacency_t contain the connection information between components
    ///
    /// adjacency_t is included by molecule_t for the recording of connections
    /// between components. As a matter of fact, objects in a component are 
    /// identified by their ID numbers, thus adjacency_t also use this ID number 
    /// to identify connection.
    ///
    /// \sa component_t, molecule_t
    class mcmprela_t
    {
      public:

        /// \brief iterator used by adjacency, which is acturally 
        ///        vector< int >::iterator
        typedef vector< int >::iterator iterator;

        /// \brief construction function
        ///
        mcmprela_t( );

        /// \brief copy construction function 
        mcmprela_t( mcmprela_t const& rhs );
        
        /// \brief deconstruction function
        ~mcmprela_t( );
        
        void swap(mcmprela_t& rhs);

        /// \brief make connection between objects
        /// \param a  ID number of an object
        /// \param b  ID number of an object
        /// \return true if connection made successfully, 
        ///         false if a and b are already connected.
        bool add( int a, int b );

        /// \brief test if connection exist between a and b
        /// \param a  ID number of an object
        /// \param b  ID number of an object
        /// \return test result if there is a connection between a and b.
        bool has( int a, int b ) const;
        
        /// \brief remove connection between a and b
        /// \param a  ID number of an object
        /// \param b  ID number of an object
        /// \return true if connection removed successfully, 
        ///         false if there is no connection between a and b.
        bool remove( int a, int b );

        /// remove all connections of a
        void clear( int a );

        /// \brief get iterator pointing to the connection.
        /// \param a  ID number of an object
        /// \param b  ID number of an object
        /// \return the iterator pointing to the connection, usually this 
        ///         function is used with begin() to give the ID the connection.
        /// \sa     begin()
        iterator find ( int a, int b ) const;

        /// \brief starting point of an object's connection
        /// \param id ID number of an object
        iterator begin( int id ) const;
        
        /// \brief ending point of an object's connection
        /// \param id ID number of an object
        iterator end  ( int id ) const;

        /// \brief number of an object's connection
        int size( int id ) const;

        int size( void ) const 
        {
            return m_content.size();
        }

        bool empty( void ) const;

        vector<int>& getvec( int a );

        const vector<int>& getvec( int a ) const;
        
      private:
        
        vector< vector<int> > m_content;

    }; // end class adjacency_t

} // namespace mort

#endif

#ifndef MORTSRC_OBJECT_COMP_HPP
#define MORTSRC_OBJECT_COMP_HPP

#include <map>
#include <boost/shared_ptr.hpp>
#include "root.hpp"
#include "atmv.hpp"

namespace mort
{
    using std::pair;
    using boost::shared_ptr;

    /// \brief database_t store grouped entity_t together
    /// \ingroup objdef
    /// database_t provides a mechanism to group entities together, and to visit
    /// them by name. Entities could be any type inherited from entity_t, including 
    /// subset_t, molecule_t, and database_t. Objects such atom, bond, resd, etc 
    /// can not by put into database because they are not entities, they do not have 
    /// their own storage space.
    /// \sa entity_t, subset_t, molecule_t 
    ///
    //
    class morf_t;
    class atmvec;
    class namemap_t;
    class molecule_t;
    class database_t;

    typedef shared_ptr<root_t> entity_ptr;
    typedef shared_ptr<namemap_t> namemap_ptr;
    typedef shared_ptr<molecule_t> molecule_ptr;
    typedef shared_ptr<database_t> database_ptr;

    class database_t : public root_t
    {
    public:

        /// \brief iterator type of database
        typedef vector< pair<string, entity_ptr> >::iterator iterator;

        /// \brief const iterator type of database
        typedef vector< pair<string, entity_ptr> >::const_iterator const_iterator;

        /// \brief constructor
        database_t(); 
       
        virtual ~database_t();
       
        /// \brief query an entity by name
        bool has( const string& name ) const;

        /// \brief get the entity with given name
        entity_ptr get( const string& name ) const;

        molecule_ptr get_mol(const string& name) const;

        database_ptr get_mdb(const string& name) const;

        atmvec_ptr get_avec(const string& name) const;

        bndvec_ptr get_bvec(const string& name) const;

        namemap_ptr get_nmap(const string& name) const;
 
	morf_t get_atom(const string& mask) const;

        morf_t get_resd(const string& mask) const;
        
        /// \brief set the entry name to an entity
        void set(const string& name, const entity_ptr& pmol );

        /// \brief set the entry name to an entity
        void add(const string& name, const entity_ptr& pmol );

        /// \brief remove an entry with given name.
        bool remove(const string& name);

        /// \brief the begin iterator of entries.
        iterator begin()
        {
            return m_data.begin();
        }

        /// \brief the ending iterator of entries
        iterator end()
        {
            return m_data.end();
        }

        /// \brief the begin constant iterator of entries
        const_iterator begin() const
        {
            return m_data.begin();
        }

        /// \brief the ending constant iterator of entries
        const_iterator end() const
        {
            return m_data.end();
        }

        int size() const
        {
            return m_data.size();
        }

    private:

        vector< pair<string, entity_ptr> > m_data;
    };
    
    typedef shared_ptr<database_t> database_ptr;

    database_t& mortenv();

} // namespace mort

#endif


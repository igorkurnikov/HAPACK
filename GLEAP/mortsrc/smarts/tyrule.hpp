#ifndef MORTSRC_SMARTS_TYRULE_HPP
#define MORTSRC_SMARTS_TYRULE_HPP

#include <map>
#include <vector>
#include <string>
#include <object.hpp>
#include "subst.hpp"

namespace mort
{

    struct item_t : public std::vector<std::string>
    {
        item_t( int size );
        
        void read_pattern( int id );

        molecule_t pattern;
    };
    
    class typing_rule : public std::vector<item_t>
    {
    public:
        
        typing_rule( const char* filename );
    
        virtual ~typing_rule();

        void read_title( std::istream& stream );
    
        void read_item( std::istream& stream );
        
        const char* get(const morf_t& atom, const char* name) const;
     
        const char* get(int i, const char* name) const;
        
    private:

        std::map<int,int> m_title;

    };

    bool read_smiles( const char* exp, molecule_t& mol );

    bool read_smarts( const char* exp, molecule_t& mol );


 
    
} // namespace mort


#endif

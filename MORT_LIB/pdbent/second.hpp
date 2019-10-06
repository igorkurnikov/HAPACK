#ifndef MORTSRC_PDBENT_SECOND_HPP
#define MORTSRC_PDBENT_SECOND_HPP

namespace mort
{
    class molecule_t;    

    namespace pdbent
    {
        
        struct second_s
        {
            second_s()
            {
                start = " ";
                end   = " ";
            }

            string start;
            string end;
            int type;
        };
 
        void read_second( istream& is, second_s& second );
 
        void apply_second( const second_s& second, molecule_t& mol );
       
    } // namespace pdbent

} // namespace mort

#endif

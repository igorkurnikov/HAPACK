#include <iostream>
#include <object.hpp>
#include <guilib.hpp>

namespace mort
{
    using std::find_if;

    namespace pdbent
    {

	void read_bond( istream& is, molecule_t& m )
	{
            const int SEQ_START = 6;
            const int SEQ_WIDTH = 5;
	    string line;
            std::getline( is, line );

            int maxseq = std::min( static_cast<size_t>(4), (line.length()-SEQ_START)/SEQ_WIDTH );
            if( maxseq <= 0 )
            {
                std::cout << "Warning: cannot recognize the following PDB CONECT card: " << std::endl;
                std::cout << "Warning:   "  << line << std::endl;
                std::cout << "Warning: it will be ignored" << std::endl;
                return;
            }

            vector<int> seqs;
            for( int i=0; i < maxseq; ++i )
            {
                int s = atoi( line.substr(SEQ_START+SEQ_WIDTH*i, SEQ_WIDTH).c_str() );
                if( s!=0 ) seqs.push_back(s);
            }

	    atom_t a1 = *find_if( m.atom_begin(), m.atom_end(), iparm_cmper1(SEQUENCE, seqs[0]) );
            for( unsigned int i=1; i < seqs.size(); ++i )
            {
		atom_t a2 = *find_if( m.atom_begin(), m.atom_end(), iparm_cmper1(SEQUENCE, seqs[i]) );
		if( bond_t::has(a1, a2) )
		{
		    continue;
		}
                    
                bond_t b = bond_t::create( a1, a2 );
                b.set_i( ORDER, 1 );
            }
        }

        //
        bool is_regular_resd( const string& name )
        {
            static char NAMES[][10] = { "ACE", "ALA", "ARG","ASH","ASN","ASP","CYM","CYS","CYX","GLH","GLN","GLU","GLY",
            "HID","HIE","HIP","ILE","LEU","LYN","LYS","MET","NHE","NME","PHE","PRO","SER","THR","TRP","TYR","VAL","HIS", 
	    "CALA","CARG","CASN","CASP","CCYS","CCYX","CGLN","CGLU","CGLY","CHID","CHIE","CHIP","CILE","CLEU","CLYS","CHIS",
	    "CMET","CPHE","CPRO","CSER","CTHR","CTRP","CTYR","CVAL","NALA","NARG","NASN","NASP","NCYS","NCYX","NGLN","NHIS",
            "NGLU","NGLY","NHID","NHIE","NHIP","NILE","NLEU","NLYS","NMET","NPHE","NPRO","NSER","NTHR","NTRP","NTYR",
	    "NVAL", "WAT","HOH","DA","DA3","DA5","DAN","DC","DC3","DC5","DCN","DG","DG3","DG5","DGN","DT","DT3","DT5","DTN",
            "RA","RA3","RA5","RAN","RC","RC3","RC5","RCN","RG","RG3","RG5","RGN","RU","RU3","RU5","RUN"};

	    static int NNAME = sizeof(NAMES)/10/sizeof(char);

	    for( int i=0; i < NNAME; ++i )
	    {
	        if( name==NAMES[i] )
	            return true;
	    }

	    return false;
        }


        //
        //  Two types of bonds need to be marked up: 
        //    (1) inter-residue bonds
        //    (2) bond in irregular residues, or inhibitors
        //
        void mark_bond( const molecule_t& m, set<int>& conns )
        {
            funstack_t::push( "pdbent::mark_bond" );

            bonditer_t bi = m.bond_begin();
            bonditer_t be = m.bond_end();
            for( ; bi != be; ++bi )
            {
                atom_t a1 = atom_1st( *bi );
                atom_t a2 = atom_2nd( *bi );
                int id1 = a1.get_i(ID);
	        int id2 = a2.get_i(ID);

                string restype;
	        if( m.nresd() > 0 )
	        {
	            resd_t r1 = a1.resd();
	            resd_t r2 = a2.resd();

                    int rid1 = r1.get_i(ID);
                    int rid2 = r2.get_i(ID);
                    int diff = std::abs(rid1-rid2);
                    if( diff != 0 && diff != 1 )
                    {
	    	        conns.insert(id1);
		        conns.insert(id2);
                        continue;
                    }

	            restype = r1.get_s(TYPE);
                }
	        else
	        {
	            restype = m.get_s(TYPE);
	        }

                if( !is_regular_resd(restype) )
                {
		    conns.insert(id1);
		    conns.insert(id2);
                }
            }

            funstack_t::pop();
        }

    } // namespace pdbent

    void disulfide( molecule_t& m, double cutoff, console_t* pcon )
    {
        vector< atom_t > sgs;
        vector< resd_t > cys;

        resditer_t ri = m.resd_begin();
        resditer_t re = m.resd_end();
        for( ; ri != re; ++ri )
        {
            string type = ri->get_s( TYPE );
            
            if( type=="CYS" || type=="CYX" )
            {        
                sgs.push_back( atom_t::get(*ri, string("SG")) );
                cys.push_back( *ri );
            }
        }

        for( unsigned int i=0; i < sgs.size(); ++i )
        {
            for( unsigned int j=i+1; j < sgs.size(); ++j )
            { 
                if( bond_t::has( sgs[i], sgs[j] ) )
                {
                    cys[i].set_s( TYPE, "CYX" );
                    cys[j].set_s( TYPE, "CYX" );
                    break;
                }

                if( dist(sgs[i], sgs[j]) > cutoff )
                {
                    continue;
                }

                bool fix = false;

                if( pcon == NULL )
                {
                    fix = true;
                }
                else
                {
                    pcon->print( "create disulf bond between " + cys[i].get_s(NAME) + " " + cys[j].get_s(NAME) + "?(y/n) " );
 
                    if( pcon->getchar() == 'y' )
                    {
                        fix = true;
                    }
                } 

                if( fix )
                {
                    cys[i].set_s( TYPE, "CYX" );
                    cys[j].set_s( TYPE, "CYX" );
                    bond_t::create( sgs[i], sgs[j] ).set_i(ORDER, 1);
                    std::cout << "make bond between atom " << sgs[i].absid() << " " << sgs[j].absid() << std::endl;
                }
            }
        }
    } 

} // namespace mort



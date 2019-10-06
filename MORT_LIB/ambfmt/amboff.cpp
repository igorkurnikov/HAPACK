
#include <map>
#include <sstream>
#include <common.hpp>
#include <object.hpp>

namespace mort
{

    namespace amboff
    {
        class reader_i
        {
        public:
            reader_i() {};

	    virtual ~reader_i() {};

            virtual void run( istream& stream, molecule_t& mol ) const = 0;
        };

        class reader_set
        {
        private:

            reader_set() {}
            
            virtual ~reader_set() {}
            
        public:

            static void insert( const string& type, const reader_i* preader )
            {
                content().insert( make_pair(type, preader) ).second;            
            }
            
            static void read( const string& type, std::istream& stream, molecule_t& mol )
            {
                std::map< string, const reader_i* >::iterator i = content().find( type );
            
                if( i == content().end() )
                {
                    throw std::runtime_error( "unknown tag " + string(type) + "!" );
                }
                
                i->second->run( stream, mol );
            }


        private:

            static std::map< string, const reader_i* >& content()
            {
                static std::map< string, const reader_i* > inst;
                return inst;
            }
        };

        class null_reader : public reader_i
        {
        public:
            
            null_reader( )
            {
                reader_set::insert( "atomspertinfo", this );
                reader_set::insert( "childsequence", this );
                reader_set::insert( "hierarchy", this );
                reader_set::insert( "name", this );
                reader_set::insert( "residuesPdbSequenceNumber", this );
                reader_set::insert( "velocities", this );
            }

            virtual void run( istream& stream, molecule_t& ) const
            {
                while( stream && stream.peek() != '!' )
                {
                    stream.ignore(MAX_LINE_WIDTH, '\n');
                }
            }   
        };

        class resd_reader : public reader_i
        {
        public:

            resd_reader()
            {
                reader_set::insert( "residues", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
                int start = 0;
                int rprev = -1;
                while( stream && stream.peek() != '!' )
                {
                    string name;
                    int seq, childseq, end;
                    
                    stream >> name >> seq >> childseq  >> end;
                    stream.ignore(MAX_LINE_WIDTH, '\n' );

                    for( int i = start; i < end - 1; ++i )
                    {
                        assert( i >= 0 && rprev >= 0 ); 
                        atom_t atom( mol, i );
                        resd_t resd( mol, rprev );
                        resd.connect(atom);
                        atom.connect(resd);
                    }
                    
                    start = end - 1;
                    resd_t curt( mol, rprev+1 );
                    string type = name.substr( 1, name.length() - 2 );
                    curt.set_s(NAME, type);
                    curt.set_s(TYPE, type);
                    rprev += 1;
                }
                  
                for( int i = start; i < mol.natom(); ++i )
                {
                    assert( i >= 0 && rprev >= 0 ); 
                    atom_t atom(mol, i);
                    resd_t resd(mol, rprev);
                    resd.connect(atom);
                    atom.connect(resd);
                }
            }
        };

        class box_reader : public reader_i
        {
        public:
            
            box_reader()
            {
                reader_set::insert( "boundbox", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
                double flag, theta, x, y, z;
                
                stream >> flag >> theta >> x >> y >> z;
                
                if( flag > 0 )
                {
                    if( theta <= M_PI )
                    {
                        theta = theta/M_PI*180.0;
                    }

                    mol.set_v(BOX, makevec(x, y, z, theta) );
                    mol.set_i(SOLUTE, BOX);
                }
                
                stream.ignore(MAX_LINE_WIDTH, '\n' );
                assert( stream.peek() == '!' );
            }
        };    

        class cap_reader : public reader_i
        {
        public:
            
            cap_reader()
            {
                reader_set::insert( "solventcap", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
                double flag, x, y, z, radius;
                
                stream >> flag >> x >> y >> z >> radius;
                
                if( flag > 0 )
                {
                    mol.set_i(SOLUTE, CAP);
                    mol.set_v(CAP, makevec(x, y, z, radius) );
                }
                
                stream.ignore( 160, '\n' );
                assert( stream.peek() == '!' );
            }
        };    

        class atom_reader : public reader_i
        {
        public:

            atom_reader( )
            {
                reader_set::insert( "atoms", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const 
            {
                while( stream.peek() != '!' )
                {
                    string line;
                    std::getline( stream, line );

                    int    typex, resx, flags, seq, element;
                    string name, type, poltype;
                    double pchg;
                    atom_t atom = mol.create(ATOM);                    
                    
                    std::istringstream ls( line );
                    ls >> name >> type >> typex >> resx >> flags >> seq >> element >> pchg >> poltype;

                    atom.set_s(NAME, strip_quota(name) );
                    atom.set_s(TYPE, strip_quota(type) );
                    atom.set_d(PCHG, pchg );
                    atom.set_i(ELEMENT, element);

                    if( !empty(poltype) )
                    {
                        atom.set_s(POLTYPE, poltype);
                    }
                }
            }
        };
            
        class bond_reader : public reader_i
        {
        public:

            bond_reader( )
            {
                reader_set::insert( "connectivity", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const 
            {
                while( stream && stream.peek() != '!' )
                {
                    int first, second, order;
		    stream >> first >> second >> order;
		    stream.ignore(MAX_LINE_WIDTH, '\n' );

		    atom_t a1(mol, first-1);
		    atom_t a2(mol, second-1);
                    bond_t b = bond_t::create(a1,a2);
                    b.set_i(ORDER,order);
                }
            }
        };

        class unitconn_reader : public reader_i
        {
        public:

            unitconn_reader( )
            {
                reader_set::insert( "connect", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
                int head, tail;
                
                stream >> head;
                stream >> tail;
                stream.ignore(MAX_LINE_WIDTH, '\n' );
                    
                mol.set_i(HEAD, head);
                mol.set_i(TAIL, tail);

                assert( stream.peek() == '!' );
            }
        };

        class resdconn_reader : public reader_i
        {
        public:

            resdconn_reader( )
            {
                reader_set::insert( "residueconnect", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
                while( stream && stream.peek() != '!' )
                {
                    int head, tail;                
                    stream >> head;
                    stream >> tail;
                    stream.ignore(MAX_LINE_WIDTH, '\n' );

                    resd_t resd = resd_t::create( mol );                    
                    resd.set_i(HEAD, head);
                    resd.set_i(TAIL, tail);
                }

                assert( stream.peek() == '!' );
            }
        };
        
        class position_reader : public reader_i
        {
        public:

            position_reader( )
            {
                reader_set::insert( "positions", this );
            }
            
            virtual void run( istream& stream, molecule_t& mol ) const
            {
	        atomiter_t atom = mol.atom_begin();

		for( ; atom != mol.atom_end(); ++atom )
		{
		    read_vparm(stream, *atom, POSITION, 3);
		    stream.ignore(MAX_LINE_WIDTH, '\n' );
		}

                assert( stream.peek() == '!' );
            }
        };

        bool read_head( istream& stream )
        {
            std::ios::pos_type pos = stream.tellg( );

            string title;
            
            getline( stream, title );

            if( title != "!!index array str" )
            {
                stream.seekg( pos );
                return false;
            }
            
            while( stream && stream.peek() != '!' )
            {
                getline( stream, title );
            }

            return true;
        }

        bool read_section( istream& stream, molecule_t& mol )
        {
            std::ios::pos_type pos = stream.tellg( );

            string entry = next_word( stream, '.' );
            stream.ignore();
            
            assert( entry == "!entry" );

            string molname = next_word( stream, '.' );
            stream.ignore();

			string curname;
			if( mol.get_s(NAME, curname) ) 
			{
				if(curname != molname)
				{
					stream.seekg( pos );
					return false;
				}
			}
			else
			{
				mol.set_s(NAME, molname);
				mol.set_s(TYPE, molname);
			}
                        
            string unit = next_word( stream, '.' );
            stream.ignore();
            
            assert( unit == "unit" );

            string info = next_word( stream );
            stream.ignore(MAX_LINE_WIDTH, '\n' );
            
            reader_set::read( info, stream, mol );
            return true;
        }

        void write_atom( ostream* os, const morf_t& atom )
        {
            *os << " ";
            *os << "\"" << atom.get_s(NAME) << "\" ";
            *os << "\"" << atom.get_s(TYPE) << "\" ";
            *os << "0 ";  // typex
            *os << atom.get_i(RESID) << " ";
            *os << "131072 "; // flags

            if( atom.nresd() > 0 )
            {
                *os << atom.absid() - atom.resd_begin()->atom_begin()->absid() + 1 << " ";
            }
            else
            {
                *os << atom.absid() + 1 << " ";
            }

            *os << atom.get_i(ELEMENT) << " ";
            *os << format("%.6f") % atom.get_d(PCHG) << " "; 

	    *os << std::endl;
        }
            

        void write_atoms( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.atoms table  ";
            os << "str name  str type  int typex  int resx  ";
            os << "int flags  int seq  int elmnt dbl chg str" << std::endl;

            atomiter_t ai = mol.atom_begin();
            atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                write_atom( &os, *ai );
            }
        }
        
        void write_boundbox( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.boundbox array dbl" << std::endl;
            
	    numvec box(4);
            if( mol.get_v(BOX, box) )
            {
                os << " 1.000000" << std::endl;
                os << " " << box[3]/180.0*M_PI << std::endl;
                os << " " << box[0] << std::endl;
                os << " " << box[1] << std::endl;
                os << " " << box[2] << std::endl;
            }
            else
            {
                os << " -1.000000" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
            }
        }
 
        void write_childseq( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.childsequence single int" << std::endl;
            
            os << " " << mol.nresd() + 1 << std::endl;
        }
        
        void write_connect( ostream& os, const molecule_t& mol )
        {
	    int head,tail;
            os << "!entry." << mol.get_s(NAME) << ".unit.connect array int" << std::endl;
            os << " " << ( mol.get_i(HEAD,head) ? head : 0 ) << std::endl;
            os << " " << ( mol.get_i(TAIL,tail) ? tail : 0 ) << std::endl;
        }

        void write_bonds( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;

            atomiter_t ai = mol.atom_begin();
            atomiter_t ie = mol.atom_end();
            for( ; ai != ie; ++ai )
            {
                atomiter_t aj = ai->atom_begin();
                atomiter_t je = ai->atom_end();
                for( ; aj != je; ++aj )
                {
                    int id1 = ai->get_i(ID);
                    int id2 = aj->get_i(ID);
                    if( id1 < id2 )
                    {
                        os << " " << ai->get_i(ID);
                        os << " " << aj->get_i(ID);
                        os << " " << 1 << std::endl;
                    }
                }
            }            
        }
        
        void write_hierarchy( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << std::endl;   

            resditer_t resd = mol.resd_begin();
            for( ; resd != mol.resd_end(); ++resd )
            {
                os << " ";
                os << "\"U\" 0 ";
                os << "\"R\" " << resd->get_i(ID) << std::endl;
                
                atomiter_t atom = resd->atom_begin();
                for( ; atom != resd->atom_end(); ++atom )
                {
                    os << " ";
                    os << "\"R\" " << resd->get_i(ID) << " ";
                    os << "\"A\" " << atom->get_i(ID) << std::endl;
                }
            }
        }

        void write_unit_name( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.name single str" << std::endl;
            
            os << " \"" << mol.get_s(NAME) << "\"" << std::endl;
        }
        
        void write_position( ostream* os, const morf_t& atom )
        {
            numvec pos = atom.get_v(POSITION);
            *os << " ";
            *os << format( " %10.6f" ) % pos[0] << " ";
            *os << format( " %10.6f" ) % pos[1] << " ";
            *os << format( " %10.6f" ) % pos[2] << std::endl;
        }

        void write_pertinfo( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;

            atomiter_t ai = mol.atom_begin();
            atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                os << " \"" << ai->get_s(NAME) << "\"";
                os << " \"" << ai->get_s(TYPE) << "\"";
                os << " 0 -1 0.0" << std::endl;
            }
        }

        void write_positions( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
            atomiter_t ai = mol.atom_begin();
            atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                write_position( &os, *ai );
            }
        }
        
        void write_resdconn( ostream* os, const morf_t& resd )
        {
	    int head, tail, conn;
            *os << " ";
            *os << ( resd.get_i(HEAD,head) ? head : 0 ) << " ";
            *os << ( resd.get_i(TAIL,tail) ? tail : 0 ) << " ";
            *os << ( resd.get_i(CONN1,conn) ? conn : 0 ) << " ";
            *os << "0 0 0" << std::endl;
        }

        void write_resdconns( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
            
            resditer_t ri = mol.resd_begin();
            resditer_t re = mol.resd_end();
            for( ; ri != re; ++ri )
            {
                write_resdconn( &os, *ri );
            }
        }
        
        void write_resd( ostream* os, const morf_t& resd )
        {
            *os << " \"" << resd.get_s(TYPE) << "\" ";
            *os << resd.get_i(ID) << " ";
            *os << resd.natom() + 1 << " ";
            *os << resd.atom_begin()->get_i(ID) << " ";
            *os << "\"p\" "; // resd type p for protein
            *os << 0 << std::endl; // imaging atom, do not know what does it mean
        }

        void write_resds( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.residues table  ";
            os << "str name  int seq  int childseq  int startatomx  str restype int imagingx" << std::endl;
            
            resditer_t ri = mol.resd_begin();
            resditer_t re = mol.resd_end();
            for( ; ri != re; ++ri )
            {
                write_resd( &os, *ri );
            }
        }
        
        void write_pdbseq( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.residuesPdbSequenceNumber array int" << std::endl;
            for( int i=0; i < mol.nresd(); ++i )
            {
                os << " " << i+1 << std::endl;
            }
        }
        
        void write_solventcap( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.solventcap array dbl" << std::endl;
            
	    numvec cap(4);
            if( mol.get_v(CAP, cap) )
            {
                os << " 1.000000" << std::endl;
                os << " " << cap[0] << std::endl;
                os << " " << cap[1] << std::endl;
                os << " " << cap[2] << std::endl;
                os << " " << cap[3] << std::endl;
            }
            else
            {
                os << " -1.000000" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
                os << " 0.0" << std::endl;
            }
        }        
        
        void write_velocity( ostream& os, const molecule_t& mol )
        {
            os << "!entry." << mol.get_s(NAME) << ".unit.velocities table  ";
            os << "dbl x  dbl y  dbl z" << std::endl;

            for( int i=0; i<mol.natom(); ++i )
            {
                os << " 0.0 0.0 0.0" << std::endl;
            }
        }

    } // namespace amber_off

    bool read_off( istream& stream, molecule_t& mol )
    {
        molecule_t temp;
        
        amboff::read_head( stream );
 
        while( stream.peek() == '!' && amboff::read_section( stream, temp ) );

        mol.swap( temp );
        
        return true;
    }

    void write_off( ostream& os, const molecule_t& mol )
    {
        amboff::write_atoms( os, mol );

        amboff::write_pertinfo( os, mol );

        amboff::write_boundbox( os, mol );
        
        amboff::write_childseq( os, mol );
        
        amboff::write_connect( os, mol );
        
        amboff::write_bonds( os, mol );
        
        amboff::write_hierarchy( os, mol );
        
        amboff::write_unit_name( os, mol );
        
        amboff::write_positions( os, mol );
        
        amboff::write_resdconns( os, mol );
        
        amboff::write_resds( os, mol );
        
        amboff::write_pdbseq( os, mol );
        
        amboff::write_solventcap( os, mol );
        
        amboff::write_velocity( os, mol );
    }


} // namespace mort

mort::amboff::null_reader g_null_reader;

mort::amboff::atom_reader g_atom_reader;

mort::amboff::bond_reader g_bond_reader;

mort::amboff::resd_reader g_resd_reader;

mort::amboff::box_reader g_box_reader;

mort::amboff::cap_reader g_cap_reader;

mort::amboff::position_reader g_position_reader;

mort::amboff::resdconn_reader g_resdconn_reader;

mort::amboff::unitconn_reader g_unitconn_reader;



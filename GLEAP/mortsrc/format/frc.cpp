#include <sstream>
#include <stdexcept>
#include <common.hpp>
#include <object.hpp>

namespace mort
{
    using namespace boost;

    using std::map;
            
    using std::istringstream;
      
    using std::runtime_error;

    typedef void (*sectfunc_t)( atmvec&, const vector<string>& );

    string str( const morf_t& mo );
 
    class frcmod_counter
    {
    public:

        static int& get()
        {
            static int nfrcmod;
            return nfrcmod;
        }

        static void init( molecule_t& ff )
        {
            if( ff.has_i(FRCMOD) )
            {
                get() = ff.get_i(FRCMOD)+1;
            }
            else
            {
                get() = 0;
            }

            ff.set_i(FRCMOD, get());
        }
    };

    void read_sect(istream& is, molecule_t& ff, int natom, sectfunc_t sectfunc, const string& type)
    {
        atmvec as(natom, atom_t(ff,-1));
        string line;
        while( getline(is, line) && !empty(line) )
        {
            if( line[0]=='#' )
            {
                continue;
            }

            if( type!="ATOM" && type!="VDWP" ) replace( line, '-', ' ' );

            vector<string> ts = split( line, " \t" );
            assert( ts.size()>0 && ts[0].length()>0 );

            if( !isdigit(ts[0][0]) )
            {
                if( (int)ts.size() < natom )
                { 
                    throw runtime_error( "Error: no enough atoms were specified for " + type + "\nInput: " + line );
                }

                bool ignore = false;
                for( int i=0; i < natom; ++i )
                {
                    atom_t a(ff, -1);
                    if( atom_t::get( ff, ts[i], a) )
                    {
                        as[i] = a;
                    }
                    else if( type=="ATOM" )
                    {
                        as[i] = atom_t::create( ff, ts[i] );
                    }
                    else
                    {
                        ignore = true;
                        std::cout << "Info: ignoring unknown atom type: " << ts[i] << std::endl;
                        std::cout << "Info: while reading " << type <<  ". Input: " << line << std::endl;
                    }
                }

                if( ignore ) 
                {
                    continue;
                }

                if( type=="DIHE" )
                {
                    as.set_d(SCEE, 1.2);
                    as.set_d(SCNB, 2.0);
                    for( unsigned int i=0; i < ts.size(); ++i )
                    {
                        if( ts[i].length() <= 5 ) continue;
            
                        string tit = ts[i].substr(0, 5);
                        if( tit=="SCEE=" )
                        {
                            double scee = atof( ts[i].c_str()+5 );
                            as.set_d( SCEE, scee );
                            //std::cout << "setting scee to: " << scee << std::endl;
                        }

                        if( tit=="SCNB=" )
                        {
                            double scnb = atof( ts[i].c_str()+5 );
                            as.set_d( SCNB, scnb );
                            //std::cout << "setting scnb to: " << scnb << std::endl;
                        }
                    }
                }
            }

            sectfunc( as, ts );
        }

    }


    void skip_hydr( istream& is )
    {
        string line;
        getline( is, line ); // && count_item( line ) >= 20 );
    }

    void skip_hbnd( istream& is )
    {
        string line;
        while( getline(is, line) && !empty(line) )
		;
    }

    void read_atom( atmvec& as, const vector<string>& ts )
    {
        assert( as.size()==1 && ts.size() >= 2 );
      
        double mass = atof( ts[1].c_str() );
        double plar = 0.0;
        if( ts.size()>2 && isdigit(ts[2][0]) )
        {
            plar = atof( ts[2].c_str() );
        }

        as[0].set_d(MASS,  mass);
        as[0].set_d(POLAR, plar);
    }        

   
    void read_bond( atmvec& as, const vector<string>& ts )
    {
        assert( as.size()==2 && ts.size()>=4 );

        double force = atof( ts[2].c_str() );
        double equil = atof( ts[3].c_str() );
       
        bond_t b = bond_t::frcget( as[0], as[1] );
        b.set_d(EQUIL, equil);
        b.set_d(FORCE, force);
    }

    void read_angl( atmvec& as, const vector<string>& ts )
    {
        assert( as.size()==3 && ts.size()>=5 );

        double force = atof( ts[3].c_str() );
        double equil = atof( ts[4].c_str() );

        angl_t an = angl_t::frcget( as[0], as[1], as[2] );
        an.set_d(EQUIL, float(equil) * DEG2RAD );
        an.set_d(FORCE, force );
    }

    bool verbose()
    {
        static string dummy;
        static bool verbose = mortenv().get_s("ffverbose", dummy) && dummy=="on";
        return verbose;
    }

    void read_dihe( atmvec& as, const vector<string>& ts )
    {
        int nfrcmod = frcmod_counter::get();

        assert( ts.size() > 0  );
        double force;
        double equil;
        double period;
        double divide = 1.0;
        if( isdigit(ts[0][0]) )
        {
            assert( ts.size() >= 4 );
            divide = atof(ts[0].c_str());
            force  = atof(ts[1].c_str());
            equil  = atof(ts[2].c_str());       
            period = atof(ts[3].c_str());
        }
        else
        {
            assert( ts.size() >= 8 );
            divide = atof(ts[4].c_str());
            force  = atof(ts[5].c_str());
            equil  = atof(ts[6].c_str());       
            period = atof(ts[7].c_str());
        }

        if( divide==0.0 )
        {
            divide = 1.0;
        }

        int p = abs(int(floor(period+0.5)));
                    
        dihvec dihs = dihe_t::get(as[0], as[1], as[2], as[3]);

        for( unsigned int i=0; i < dihs.size(); ++i )
        {
            if( dihs[i].get_i(FRCMOD) < nfrcmod && dihs[i].get_d(FORCE) > 1e-6 )
            {
                dihs[i].set_d(FORCE, 0.0);           
                if( verbose() )
                {
                    std::cout << "NOTE: torsion term: " << str(dihs[i]) << " is overwritten by " << std::endl;
                    std::cout << "NOTE:         term: ";
                    std::cout << as[0].name() << "-" << as[1].name() << "-";
                    std::cout << as[2].name() << "-" << as[3].name() << " ";
                    std::cout << "v" << p << " in frcmod file." << std::endl;
                }
            }
        }

        bond_t::frcget( as[1], as[2] );
                    
        dihe_t di = dihe_t::frcget(as[0],as[1],as[2],as[3],p);

        di.set_d(FORCE,  force/divide);
        di.set_d(EQUIL,  equil*DEG2RAD);
        di.set_i(PERIOD, p );
        di.set_i(FRCMOD, nfrcmod );
        di.set_d(SCEE,   as.get_d(SCEE) );
        di.set_d(SCNB,   as.get_d(SCNB) );
    }

    void read_impr( atmvec& as, const vector<string>& ts )
    {
        assert( as.size()==4 && ts.size()>=7 );

        double divide, force, equil, period;
        if( ts.size()>=8 && isdigit(ts[7][0]) )
        {
            divide = atof( ts[4].c_str() );
            force  = atof( ts[5].c_str() );
            equil  = atof( ts[6].c_str() );
            period = atof( ts[7].c_str() );
        }
        else
        {
            divide = 1.0;
            force  = atof( ts[4].c_str() );
            equil  = atof( ts[5].c_str() );
            period = atof( ts[6].c_str() );
        }

        if( divide== 0.0 )
        {
            divide = 1.0;
        }


        bond_t::frcget( as[2], as[3] );
        int p = abs(int(floor(period+0.5)));
        impr_t im = impr_t::frcget( as[0],as[1],as[2],as[3],p );


        im.set_d(FORCE, force);
        im.set_d(EQUIL, equil * DEG2RAD);
        im.set_i(PERIOD,p );
        im.set_d(SCEE,  0.0);
        im.set_d(SCNB,  0.0);
    }

    void read_vdwmap( istream& is, molecule_t& ff, map< string, string >& vdwmap )
    {
        string line;

        while( getline( is, line ) && ! empty( line ) )
        {
            if( line[0]=='#' ) continue;

            istringstream ls( line );            

            string name, alias;
            ls >> name;

            while( ls >> alias )
            {
		if( atom_t::has(ff, alias) )
                {
                    vdwmap[ alias ] = name;
                }
                else
                {
                    std::cout << "Info: ignoring unknown atom type: " << alias << std::endl;
                }
            }
        }

        getline( is, line ); // ignore the line "MOD4 RE" in parm99.dat, doesn't know what does it mean.
    }
            

    void read_vdwp( atmvec& as, const vector<string>& ts )
    {
        assert( as.size()==1 && ts.size()>=3 );

        double rstar = atof( ts[1].c_str() );
        double depth = atof( ts[2].c_str() );

        as[0].set_d(RSTAR, rstar);
        as[0].set_d(DEPTH, depth);
    }
               

    void assign_vdw( molecule_t& mol, const map< string, string >& vdwmap )
    {
        map< string, string >::const_iterator i = vdwmap.begin();
        for( ; i != vdwmap.end(); ++i )
        {
            atom_t dst = atom_t::get( mol, i->first );
            atom_t src = atom_t::get( mol, i->second );
                    
            dst.set_d(DEPTH, src.get_d(DEPTH) );
            dst.set_d(RSTAR, src.get_d(RSTAR) );
        }                
    }

    void read_14nb( atmvec& as, const vector<string>& ts)
    {
        double scee = atof( ts[4].c_str() );
        double scnb = atof( ts[5].c_str() );

        dihvec dihs = dihe_t::get(as[0],as[1],as[2],as[3]);
        for( unsigned int i=0; i < dihs.size(); ++i )
        {
            dihs[i].set_d( SCEE, scee );
            dihs[i].set_d( SCNB, scnb );
        } 
    } 


    void read_full_frc( istream& is, molecule_t& ff );

    void read_frc( istream& is, molecule_t& ff )
    {
        frcmod_counter::init( ff );

        string title;
        getline( is, title );

        string keyw = peek_keyw( is );

        // skip empty lines after title
        while( is && empty(keyw) )
        {
            is.ignore( MAX_LINE_WIDTH, '\n' );
            keyw = peek_keyw(is);
        }

        while( is && ! empty( keyw ) )
        {
            if( keyw == "MASS" )
            {
                is.ignore( MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 1, read_atom, "ATOM" );
            }
            else if( keyw == "BOND" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 2, read_bond, "BOND" );
            }
            else if( keyw == "ANGL" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 3, read_angl, "ANGL" );
            }
            else if( keyw == "DIHE" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 4, read_dihe, "DIHE" );
            }
            else if( keyw == "IMPR" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 4, read_impr, "OOPS" );
            }
            else if( keyw == "NONB" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 1, read_vdwp, "VDWP" );
            }
            else if( keyw == "14NB" )
            {
                is.ignore(MAX_LINE_WIDTH, '\n' );
                read_sect( is, ff, 4, read_14nb, "14NB" );
            }
            else
            {
                read_full_frc( is, ff );
                return;
            }
            
            keyw = peek_keyw( is );
        }

    }
    
    void read_full_frc( istream& is, molecule_t& ff )
    {
        frcmod_counter::init( ff );

        map< string, string > vdwmap;

        //std::cout << "reading atom" << std::endl;
        atom_t::create( ff, "X" );
        read_sect( is, ff, 1, read_atom, "ATOM" );

        skip_hydr( is );
        
        //std::cout << "reading bond" << std::endl;
        read_sect( is, ff, 2, read_bond, "BOND" );
    
        //std::cout << "reading angl" << std::endl;
        read_sect( is, ff, 3, read_angl, "ANGL" );
        
        //std::cout << "reading dihedral" << std::endl;
        read_sect( is, ff, 4, read_dihe, "DIHE" );
        
        //std::cout << "reading improper" << std::endl;
        read_sect( is, ff, 4, read_impr, "IMPR" );

        //std::cout << "reading hbond" << std::endl;            
        skip_hbnd( is );

        //std::cout << "reading vdw map" << std::endl;
        read_vdwmap( is, ff, vdwmap );
        
        //std::cout << "reading vdw" << std::endl;
        read_sect( is, ff, 1, read_vdwp, "VDWP" );

        assign_vdw( ff, vdwmap );
        
    }



}  // namespace mort



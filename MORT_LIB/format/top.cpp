#include <fstream>
#include <ambfmt.hpp>


namespace mort
{

    void load_frc( const string& file, molecule_t& ff)
    {
        std::ifstream is( file.c_str() );

        read_frc(is, ff);

	if( ff.has_i("NFRCMOD") )
	{
            int nfrcmod = ff.get_i("NFRCMOD");

            string frcmods;
            ff.get_s("FRCMODS", frcmods);

            frcmods += file;
            frcmods += "\n";
            ff.set_i("NFRCMOD", nfrcmod+1 );
            ff.set_s("FRCMODS", frcmods);
        }
        else
        {
            ff.set_i("NFRCMOD", 0);
        }
    }

    void load_amoeba_frc( const string& file, molecule_t& atomff, molecule_t& poleff )
    {
        std::ifstream is( file.c_str() );
        read_amoeba_frc( is, atomff, poleff );
    }
 
    void save_top( molecule_t& m, const string& file, const molecule_t& ffp )
    {
        std::ofstream os( file.c_str() );
        write_amber_prmtop( os, m, ffp );
    }

} // namespace mort


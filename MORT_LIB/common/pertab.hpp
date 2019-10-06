#ifndef DDUTILS_MORT_PERTAB_H
#define DDUTILS_MORT_PERTAB_H

#include <string>

namespace mort
{
    /// pertab stands for Elements PERiodic TABle. 
    /// It is used to provide information (symbol, color, etc.) for elements
    /// the data is originally from OElib 2.0


    using std::string;

    static const int NUMBER_ELEMENTS = 110;
    static const int HYDROGEN = 1;
    static const int CARBON = 6;
    static const int NITROGEN = 7;
    static const int OXYGEN = 8;
    static const int FLUORINE = 9;
    static const int SILICON = 14;
    static const int PHOSPHORUS = 15;
    static const int SULFUR = 16;
    static const int CHLORINE = 17;
    static const int BROMINE = 35;
    static const int IODINE = 53;

    class pertab_t
    {
      private:

        pertab_t( const char* buf );
    
        static const pertab_t& instance( );
    
      public:
    
        static int get_element( const string& symbol );    

        static int get_element( char symchar );

        static double get_weight( size_t element );    

        static double get_red( int element );

        static double get_blue( int element );
        
        static double get_green( int element );

        static double get_rvdw( int element );

        static const char* get_symbol( size_t element );
    
      private:

        string m_symbols[ NUMBER_ELEMENTS ];
    
        double m_weights[ NUMBER_ELEMENTS ];

        double m_colors [ NUMBER_ELEMENTS ][3];

        double m_rvdws  [ NUMBER_ELEMENTS ];
    };



} // namespace mort

#endif

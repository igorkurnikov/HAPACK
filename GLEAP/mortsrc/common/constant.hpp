#ifndef MORTSRC_COMMON_CONST_HPP
#define MORTSRC_COMMON_CONST_HPP
#include <cmath>

#ifndef M_PI
#include <math.h>
#define M_PI 3.1415926535897932385
#endif
using std::sqrt;

namespace mort
{
    /// \defgroup common Common: general used constants, types and functions
    /// \defgroup constant Constants
    /// \ingroup common
    /// @{

    static const double BONDCUT = 1.875;
    static const double BUMPCUT = 0.8;
    static const double BONDCUT2 = BONDCUT * BONDCUT;
    static const double BUMPCUT2 = BUMPCUT * BUMPCUT;

    static const double DIELFAC = 1.0 - 1.0 / 78.5;
    static const double AVERGADRO = 6.023e23;
    static const double ML_TO_A3 = 1e24;
    static const double DEG2RAD = 0.0174533;
    static const double INVCHG2 = 18.2223 * 18.2223;
    static const double PI2 = M_PI * M_PI;
    static const double TWOPI = 2.0 * M_PI;
    static const double INVSQRTPI = 1.0 / sqrt(M_PI);
    static const double DOUBLEMAX = 1.0e200;
    static const int UNKNOWN = 0; /// default value for all integer parameter
    static const int MAX_LINE_WIDTH = 4096;

    namespace daylight
    {
        enum env_e
        {
            HYBRID, DEGREE, VALENCE, HYDROGEN, NEBR
        };

        enum level_e
        {
            SINGLE, CONN, SURPRISE, COMMA, SEMICOLON
        };

        enum bracket_e
        {
            OPEN, CLOSE
        };

    }

    enum find_policy_e
    {
        FIND_ONE, FIND_ALL
    };

    enum finder_strategy
    {
        ALLOW_MULTI_VISIT = 0x0001, NO_MULTI_VISIT = 0x0000
    };

    enum predict_result
    {
        WRONG_DIRECTION, KEEP_GOING, ENOUGH
    };

    enum modifier_e
    {
        SHIFT_MASK = 1 << 0,
        LOCK_MASK = 1 << 1,
        CONTROL_MASK = 1 << 2,
        MOD1_MASK = 1 << 3,
        MOD2_MASK = 1 << 4,
        MOD3_MASK = 1 << 5,
        MOD4_MASK = 1 << 6,
        MOD5_MASK = 1 << 7,
        BUTTON1_MASK = 1 << 8,
        BUTTON2_MASK = 1 << 9,
        BUTTON3_MASK = 1 << 10,
        BUTTON4_MASK = 1 << 11,
        BUTTON5_MASK = 1 << 12,
        RELEASE_MASK = 1 << 30,
        MODIFIER_MASK = RELEASE_MASK | 0x1fff
    };

    enum quantity_e
    {
        QUANTITY_HIGH, QUANTITY_LOW
    };

} // namespace mort


#endif


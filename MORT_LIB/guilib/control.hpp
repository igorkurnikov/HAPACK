#ifndef MOLVIEW_CONTROL_H
#define MOLVIEW_CONTROL_H

#include <map>
#include <common.hpp>

namespace mort
{
    class command_i;

    /// \brief the central processing class of leap
    /// \ingroup leaplib
    class control_t
    {
    private:

        control_t();
    
        virtual ~control_t();

    public:

        /// execute a command
        static bool run( const string& command );

        /// register a command
        /// \param name command name
        /// \param name command ptr
        /// add a entry in command dictionary
        static void insert( const string& name, command_i* ptr );

        /// the begin iterator of command dictionary
        static std::map< string, command_i* >::iterator begin();

        /// the ending iterator of command dictionary
        static std::map< string, command_i* >::iterator end();
    };

    
} // namespace mort

#endif


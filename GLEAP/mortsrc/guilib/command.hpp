#ifndef MOLVIEW_ACTION_H
#define MOLVIEW_ACTION_H
#include <boost/shared_ptr.hpp>
#include <common.hpp>

namespace mort  
{
    using boost::shared_ptr;

    /// \brief command_i is the base class of all command types
    /// \ingroup leaplib
    class command_i
    {
    public:

        /// constructor
        command_i();

        /// constructor
        ///
        /// add an entry in control's command dictionary
        command_i( const string& name );

        /// deconstructor
        virtual ~command_i();

        /// command execution
        virtual bool exec() = 0;
    
        /// undoing a command execution
        virtual void undo() = 0;
    
        /// making a new command object with given argument.
        virtual shared_ptr<command_i> clone( const vector<string>& args ) const = 0;

        static void insert( const string& name, command_i* inst );

        static command_i* find( const string& name );

    private:

        static std::map<string, command_i*>& dict();
    };

    typedef vector<string> argvec_t;
    
    typedef shared_ptr<command_i> command_ptr;

} // namespace mort

#endif

#ifndef MORT_COMMON_ALGORITHM_HPP
#define MORT_COMMON_ALGORITHM_HPP

#include <string>

namespace mort
{
    using std::string;

    /// \defgroup str_algo String algortihms
    /// \ingroup common
    /// @{

    vector<string> split( const string& line, const string& dels=" " );

    /// \brief next full word in the stream before a delimeter. 
    /// 
    /// first skip spaces, then read string until space or delimer
    /// appear
    std::string next_word( std::istream& stream, char deli='\n' );

    /// \brief next alphabetical or digit in a string.
    std::string next_alnum( const char* ptr );
    
    /// \brief next digit 
    std::string next_digit( const char* ptr );

    const char* skip_alnum( const char* ptr );

    const char* skip_alpha( const char* ptr );

    const char* skip_digit( const char* ptr );

    const char* skip_float( const char* ptr );

    std::string& strip( std::string& name );
    
    void replace( std::string& str, char a, char b );
    
    string replace_copy( const std::string& str, char a, char b );

    bool empty( const std::string& str );
    
    int count_item( const std::string& line );

    string peek_keyw( std::istream& is );

    string strip_quota(const string& str);

    /// @}

    /// \defgroup generic_algorithm Generic algorithms
    /// \ingroup common
    /// @{

    /// \brief call function with different argument, sum the result together
    template<typename _InputIterator, typename _Tp, typename _UnaryOperation>
    _Tp sum(_InputIterator __first, _InputIterator __last, _Tp __init,
	       _UnaryOperation __unary_op)
    {
        for ( ; __first != __last; ++__first)
            __init += __unary_op(*__first);
        return __init;
    }

    /// \brief call function with different argument, sum the result together
    template<typename _InputIterator, typename _Tp, typename _UnaryOperation>
    _Tp sum(_InputIterator __first, _InputIterator __last, _Tp __init)
    {
        for ( ; __first != __last; ++__first)
            __init += *__first;
        return __init;
    }

    /// \brief copy elements fit for the predicate
    template<typename _InputIterator, typename _OutputIterator,
             typename _Predicate>
    _OutputIterator
    copy_if(_InputIterator __first, _InputIterator __last,
            _OutputIterator __result, _Predicate __pred)
    {
        for ( ; __first != __last; ++__first)
        {
            if (__pred(*__first))
            {
                *__result = *__first;
                ++__result;
            }
        }
        
        return __result;
    }

    /// @}

} // namespace mort

#endif

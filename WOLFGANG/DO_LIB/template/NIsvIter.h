#ifndef NIsv_ITER
#define NIsv_ITER

#include <NIsvList.h>
#include <IterBase.h>

template<class Type>
class NIsv_list_iter : private slist_base_iter {
public:
    NIsv_list_iter( NIsv_list<Type> & slist ) : slist_base_iter( &slist ) {};
    ~NIsv_list_iter() {};

    void                reset();
    void                reset( NIsv_list<Type> & );
    Type*               operator()();
    int                 operator++();
    Type*               operator+=( int adv_amount );
    NIsv_list<Type>*    container();
    Type * current();
    void append( Type * );
};

template<class Type>
inline void NIsv_list_iter<Type>::reset()
    { slist_base_iter::base_reset(); }

template<class Type>
inline void NIsv_list_iter<Type>::reset( NIsv_list<Type> & slist )
    { slist_base_iter::base_reset_list( &slist ); }

template<class Type>
Type * NIsv_list_iter<Type>::operator()()
{
    TypeLink<Type*> * tmp;
    tmp= (TypeLink<Type*> *) slist_base_iter::operator()();
    return tmp ? tmp->data : 0;
}

template<class Type>
inline int NIsv_list_iter<Type>::operator++()
    { return (slist_base_iter::operator++() != 0 ); }

template<class Type>
Type * NIsv_list_iter<Type>::operator+=( int adv_amount )
{
    TypeLink<Type*> * tmp;
    tmp= (TypeLink<Type*> *) slist_base_iter::operator+=( adv_amount );
    return tmp ? tmp->data : 0;
}

template<class Type>
inline NIsv_list<Type> * NIsv_list_iter<Type>::container()
    { return ( NIsv_list<Type> *) curr_list; }

template<class Type>
Type * NIsv_list_iter<Type>::current()
{
    Type * ret_val = 0;
    if ( curr_item == 0)
        cerr << "Iter-List-Class-Warning : There is no current item\n";
    else {
        TypeLink<Type*> * tmp;
        tmp = (TypeLink<Type*> *)(NIsv_list<Type> *) curr_item;
        ret_val = tmp->data;
    }
    return ret_val;
}

template<class Type>
void NIsv_list_iter<Type>::append( Type * datum )
{
    slink * curr_last = base_last_hit( this, curr_item );
    ((NIsv_list<Type> *)curr_list)->insert( datum);
    base_last_unhit( this, curr_last );
}
#endif

#ifndef ISV_ITER
#define ISV_ITER

#include "Isv_List.h"
#include "IterBase.h"

template<class Type>
class Isv_list_iter : private slist_base_iter {
public:
    Isv_list_iter(Isv_list<Type> &slist ) : slist_base_iter( &slist ) {};
    ~Isv_list_iter() {};

// reset iterator to first element
    void                reset();
// reset the list on which the iterator acts
    void                reset( Isv_list<Type> &  );
// advance to next element and return pointer it
    Type*               operator()();
// advance to next element and return 0 if applied to last element, else 1
    int                 operator++();
// advance forward by given amount and return pointer to element
    Type*               operator+=( int);
// return pointer to current element
    Type*               current();
// return pointer to list on which the iterator acts
    Isv_list<Type>*     container();
// insert a new element into the list after the element the iterator points to
    void                append( Type * );
};

template<class Type>
inline void Isv_list_iter<Type>::reset()
    { slist_base_iter::base_reset(); }

template<class Type>
inline void Isv_list_iter<Type>::reset( Isv_list<Type> & slist )
    { slist_base_iter::base_reset_list( &slist ); }

template<class Type>
inline Type * Isv_list_iter<Type>::operator()()
    { return (Type *) slist_base_iter::operator()(); }

template<class Type>
inline int Isv_list_iter<Type>::operator++()
    { return (slist_base_iter::operator++() != 0 ); }

template<class Type>
inline Type * Isv_list_iter<Type>::operator+=( int adv_amount )
    { return (Type *) slist_base_iter::operator+=( adv_amount ); }

template<class Type>
inline Type * Isv_list_iter<Type>::current()
    { return (Type *) curr_item; }

template<class Type>
inline Isv_list<Type> * Isv_list_iter<Type>::container()
    { return (Isv_list<Type> *) curr_list; }

template<class Type>
void Isv_list_iter<Type>::append( Type * datum )
{
    slink * curr_last = base_last_hit( this, curr_item );
    ((Isv_list<Type> *)curr_list)->insert( datum);
    base_last_unhit( this, curr_last );
}

#endif

/* Note: You have to take special care to the memory-handling!
         
         The recommended way is to allocate without giving each
         element a variable-name e.g.
                insert( new Isv_Type( initialisation ) );
         When removing an element from the list by the get operation,
         you have to delete it explicitely after using it. (Do not forget
         to keep a pointer to the element.)
         To remove and delete single elements without using them, you
         may use the destroy operation
         Finaly you can delete the elements left in the list by the
         clearAndDestroy opertion.

         Never delete directly a pointer to an element of the list, because
         your list breaks down then!!!!!!

         In order to use one element in further lists you have to use
         non-intusive lists. Then using the destroy or clearAndDestroy
         operation causes undefined values in the non-intrusive list, but
         the non-intrusive list will not break down.
         It is up to the programmer to avoid this!!!!!

         COUNTING OF ELEMENTS: The positions of the list are counted
         starting with zero. Therefore the last element has the 
         position (entries-1).

Since most of the functions are derived from the base class ListBase.h
you will find further documentaion there.
*/

#ifndef ISV_LIST
#define ISV_LIST

#include "ListBase.h"

template <class Type> class Isv_list_iter;

template<class Type>
class Isv_list : private slist_base {
public:
// add at head of list
    void                insert(Type *); 
// add at position in list
    void                include(Type *, int = 0 );
// add at tail of list
    void                append(Type *);     
// remove element from list and delete memory
    void                destroy( int = 0 );
// return pointer to element and remove it from the list
// MEMORY IS NOT DELETED
    Type*               get( int = 0 );
// return pointer to element but do not remove it
    Type*               find( int = 0 );
// return pointer to last element
    Type*               findLast();
// clear the list without deleting the memory of the objects
    void                clear();
// clear the list, but delete the memory of the objects first
    void                clearAndDestroy();
// prove if element already exists
    ListBool            contains(Type *);
// prove if list is empty
    ListBool            isEmpty();
// return number of entries
    int                 entries();
// apply function to each element
    void                forAll( void (*)(Type *, void *), void *);

    friend class Isv_list_iter<Type>;
};

template<class Type>
inline void Isv_list<Type>::insert(Type* datum)
    { slist_base::base_insert(this, datum); }

template<class Type>
inline void Isv_list<Type>::include(Type* datum, int index)
    { slist_base::base_include(this ,datum,index); }

template<class Type>
inline void Isv_list<Type>::append(Type* datum)
    { slist_base::base_append(this, datum); }

template<class Type>
void Isv_list<Type>::destroy(int index)
{
    Type* del = (Type*) slist_base::base_get(this,index,"destroy");
    if( del == 0 ) {
       cerr << "List-Class-Warning: destroy : list is already empty\n";
    }
    else {
       delete del;
    }
}

template<class Type>
inline Type* Isv_list<Type>::get(int index)
    { return (Type*) slist_base::base_get(this,index); }

template<class Type>
inline Type* Isv_list<Type>::find(int index)
    { return (Type*) slist_base::base_find(this,index); }

template<class Type>
inline Type* Isv_list<Type>::findLast()
    { return (Type*) this->last; }

template<class Type>
inline void Isv_list<Type>::clear()
    { slist_base::base_clear(this); }

template<class Type>
void Isv_list<Type>::clearAndDestroy()
{
     Type* next_link;
     Type* rover = (Type *)this->last; 
     rover = (Type *)slist_base::base_next(this, rover, True);
     while( rover != 0 ) {
         next_link = (Type *)slist_base::base_next(this, rover, False);
         delete rover;
         rover = next_link;
     }
     slist_base::base_clear(this);
}

template<class Type>
inline ListBool Isv_list<Type>::contains(Type* datum)
    { return slist_base::base_contains(this, datum); }

template<class Type>
inline ListBool Isv_list<Type>::isEmpty()
{
       if( this->last == 0 )
         return True;
       else
         return False;
}

template<class Type>
inline int Isv_list<Type>::entries()
    { return this->entry_count; }

template<class Type>
void Isv_list<Type>::forAll(register void (*func)(Type *, void *), void *datum)
{
    Type * rover = (Type *) this->last;
// the above initialization was originally done by the following line, but
// some compiler argued that rover could be used before being set.
    rover = (Type *) slist_base::base_next(this, rover, True);
    while( rover != 0) {
        (*func)(rover, datum);
        rover = (Type *) slist_base::base_next(this, rover, False);
    }
}

#endif














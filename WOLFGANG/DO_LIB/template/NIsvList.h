/*
    Be careful with the destroy and clearAndDestroy operation
    when the non-intrusive list is derived from elements of an
    intrusive list. These operations may destroy your intrusive
    list, too. Therefore it is suggested to use only the get and
    clear operation and to leave the destroying to the intrusive
    list.
*/

#ifndef NIsv_LIST
#define NIsv_LIST

#include <ListBase.h>

template<class Type>
struct TypeLink : public slink {
        Type data;
        TypeLink(const Type& datum) : data(datum) { };
};

template<class Type> class NIsv_list_iter;

template<class Type>
class NIsv_list : private slist_base {
public:
    void                insert(Type *);
    void                include(Type *, int = 0 );
    void                append(Type *);
    void                destroy( int = 0 );
    Type*               get( int = 0 );
    Type*               find( int = 0 );
    Type*               findLast();
    void                clear();
    void                clearAndDestroy();
    ListBool            contains(Type *);
    ListBool            isEmpty();
    int                 entries();
    void                forAll( void (*)(Type *, void *), void *);

    friend class NIsv_list_iter<Type>;
};

template<class Type>
inline void NIsv_list<Type>::insert( Type * datum )
    { slist_base::base_insert(this, new TypeLink<Type*>(datum)); }

template<class Type>
inline void NIsv_list<Type>::include( Type * datum, int index )
    { slist_base::base_include(this, new TypeLink<Type*>(datum), index); }

template<class Type>
inline void NIsv_list<Type>::append( Type * datum )
    { slist_base::base_append(this, new TypeLink<Type*>(datum)); }

template<class Type>
void NIsv_list<Type>::destroy(int index)
{
    TypeLink<Type*> * del_obj;
    del_obj = (TypeLink<Type*> *) slist_base::base_get(this, index,"destroy");
    if( del_obj == 0 )
       cerr << "List-Class-Warning: destroy : list is already empty\n";
    else {
       delete del_obj->data;
       delete del_obj;
    }
}

template<class Type>
Type * NIsv_list<Type>::get(int index)
{
    TypeLink<Type*> * ret_obj;
    ret_obj = (TypeLink<Type*> *) slist_base::base_get(this, index);
    Type * ret_val = ret_obj->data;
    delete ret_obj;
    return  ret_val;
}

template<class Type>
Type * NIsv_list<Type>::find(int index)
{
    TypeLink<Type*> * ret_obj;
    ret_obj = (TypeLink<Type*> *) slist_base::base_find(this, index);
    return ret_obj->data;
}

template<class Type>
inline Type * NIsv_list<Type>::findLast()
    { return ((TypeLink<Type*> *) this->last)->data; }

template<class Type>
void NIsv_list<Type>::clear()
{
     TypeLink<Type*> * next_link;
     TypeLink<Type*> * rover = 0;
     rover = (TypeLink<Type*> *) slist_base::base_next(this, rover, True);
     while( rover != 0 ) {
       next_link = (TypeLink<Type*> *) slist_base::base_next(this, rover, False);
       delete rover;
       rover = next_link;
     }
     slist_base::base_clear(this);
}

template<class Type>
void NIsv_list<Type>::clearAndDestroy()
{
     TypeLink<Type*> * next_link;
     TypeLink<Type*> * rover = 0;
     rover = (TypeLink<Type*> *) slist_base::base_next(this, rover, True);
     while( rover != 0 ) {
       next_link = (TypeLink<Type*> *) slist_base::base_next(this, rover, False);
       delete rover->data;
       delete rover;
       rover = next_link;
     }
     slist_base::base_clear(this);
}

template<class Type>
ListBool NIsv_list<Type>::contains( Type * lookup_data )
{
    TypeLink<Type*> * rover = 0;
    rover = (TypeLink<Type*> *) slist_base::base_next(this, rover, True );
    if( rover == 0 )
        return  False;
    while( rover != 0 ) {
        if( rover->data == lookup_data )
            return  True;
        rover = (TypeLink<Type*> *) slist_base::base_next(this, rover , False );
    }
    return False;
}

template<class Type>
inline ListBool NIsv_list<Type>::isEmpty()
{
        if( this->last == 0 )
           return True;
        else
          return False;
}

template<class Type>
inline int NIsv_list<Type>::entries()
    { return this->entry_count; }

template<class Type>
void NIsv_list<Type>::forAll(register void (*func)(Type *, void *), void *datum)
{
    TypeLink<Type*> * rover = 0;
    rover = (TypeLink<Type*> *) slist_base::base_next(this, rover, True);
    while( rover != 0) {
        (*func)(rover->data, datum);
        rover = (TypeLink<Type*> *) slist_base::base_next(this, rover, False);
    }
}

#endif

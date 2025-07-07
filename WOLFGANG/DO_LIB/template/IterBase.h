#ifndef ITERBASE
#define ITERBASE

#include "ListBase.h"

class slist_base_iter {
protected:
    slist_base *  curr_list;
    slink *       curr_item;
    ListBool      at_end;

// advance through list
    static inline slink * base_advance( slist_base_iter * const, int, ListBool = True);
// change pointer to last used by the derived append function
    static slink * base_last_hit( slist_base_iter * const, slink *);
// rechange pointer to last used by the derived append function
    static void base_last_unhit( slist_base_iter * const, slink *);

    inline void base_reset();
    inline void base_reset_list( slist_base * );

public:
    slist_base_iter( slist_base * list )
                          :curr_item(0), curr_list(list), at_end(False) {};
    ~slist_base_iter() {};

    inline slink * operator()() ;
    inline slink * operator++();
    inline slink * operator+=( int );
};

// This resets only the iterator to the first element

inline void slist_base_iter::base_reset() {
    at_end = False;
    curr_item = 0;
}

//  This replaces the current list by a new one

inline void slist_base_iter::base_reset_list( slist_base * list ) {
    curr_list = list;
    base_reset();
}

inline slink * slist_base_iter::operator()()
    { return base_advance (this, 1); }

inline slink * slist_base_iter::operator++()
    { return base_advance(this, 1, False); }

inline slink * slist_base_iter::operator+=(int  incr_amount )
    { return base_advance(this, incr_amount); }

inline slink * slist_base_iter::base_advance( slist_base_iter * const iterb,
                                          int adv_amount, ListBool msg )
{
    slink * list_item = iterb->curr_item;
    if( iterb->at_end ) {
       if ( msg)
           cerr << "Iter-List-Class-Warning : End of list already reached !\n";
       return 0;
    }
    if( adv_amount < 1 ) {
       cerr << "Iter-List-Class-Warning : Advance value is less then one !\n";
       adv_amount = 1;
    }
    slink * list_last = iterb->curr_list->last;
    if( list_item == 0 ) {         // set to beginning for reset list
        list_item = list_last;
    }
    while( adv_amount-- > 0 ) {
        list_item = list_item->next;
        if( list_item == list_last ) {
            iterb->at_end = True;
            if( adv_amount > 0 ) {
                cerr << "Iter-List-Class-Warning : ";
                cerr << "End of list reached before advance was complete\n";
                return 0;
            }
        }
    }
    iterb->curr_item = list_item;
    return list_item;
}

inline slink * slist_base_iter::base_last_hit( slist_base_iter * const iterb,
                                                   slink * new_last )
{
    slink * list_last = iterb->curr_list->last;
    if( list_last == new_last ) {
        list_last = 0;
    }
    iterb->curr_list->last = new_last;
    return list_last;
}

inline void slist_base_iter::base_last_unhit( slist_base_iter * const iterb,
                                                  slink * old_last )
{
    if( old_last != 0 )
       iterb->curr_list->last = old_last;
}

#endif


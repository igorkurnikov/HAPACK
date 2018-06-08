// This defines the links and the base-class of a singly-linked-list.
// The structure and most of the variable-names are adapted from Bjarne
// Stroustrup: The C++ Programming Language, Second Edition, Addison Wesley.

#include <stdlib.h>
#include "io_incl.h"

#ifndef LISTBASE
#define LISTBASE
enum ListBool {False,True};

struct slink {
    slink * next;
    slink() { next=0; }
    slink(slink* p) { next = p; }
};

class slist_base {
protected:
    slink* last;                      // last->next is head of list
    int    entry_count;               // additional to Stroustrup

// add at head of list
   static inline void base_insert(slist_base * const, slink* );
// add at abitrary position in list
   static inline void base_include(slist_base * const, slink* ,int);
// add at tail of list
   static inline void base_append(slist_base * const, slink* );
// used by base_get and base_find
   static inline int base_index_check(slist_base * const, int ,char*);
// return pointer to list element and remove the element
   static inline slink* base_get(slist_base * const, int,char*);
// basic find procedure used by base_get and base_find
   static inline slink* base_bfind(slist_base * const, int);
// return pointer to list element and do not remove
   static inline slink* base_find(slist_base * const, int);
// clear list without deleting the objects
   static inline void base_clear(slist_base * const );
//move to next element or reset to first element
   static inline slink* base_next(slist_base * const, slink* ,ListBool);
// proves if element already exists
   static ListBool base_contains(slist_base * const, slink *);

    slist_base() { last = 0; entry_count = 0; };
    slist_base(slink* datum) { last = datum->next = datum; entry_count = 1; };
    ~slist_base() {};

    friend class slist_base_iter;

};
// This defines the single list base append operation, which adds to the
// tail of the list.

inline void slist_base::base_append(slist_base * const slist, slink* datum)
{
    if (slist->last) {            // means if last already exists: (last <> 0)
       datum->next = slist->last->next;
       slist->last = slist->last->next = datum;
    }
    else
       { slist->last = datum->next = datum; }
    slist->entry_count++;
}



// This defines the single list base insert operation, which adds to the
// head of the list.

inline void slist_base::base_insert(slist_base * const slist, slink* datum)
{
    if (slist->last)        // means if last already exists: (last <> 0)
       datum->next = slist->last->next;
    else
       slist->last = datum;
    slist->last->next = datum;
    slist->entry_count++;
}

// This defines the single list base include operation, which adds at the
// specified position of the list. Note that the position is counted
// starting at zero.

inline void slist_base::base_include(slist_base * const slist, slink* datum,
                                                              int index )
{
  if( index == 0)
    { base_insert(slist,datum); }
  else if (index == slist->entry_count)
    { base_append(slist,datum); }
  else {
    base_index_check(slist,index,"include");
    slink* prev  = base_bfind(slist, index - 1);
    datum->next = prev->next;
    prev->next = datum;
    slist->entry_count++;
  }
}
// This defines the list base index check, which checks:
//               0 <= index <= entry_count-1.
// If not an error message is printed and the program is suspended

inline int slist_base::base_index_check(slist_base * const slist, int index,
                                                            char* msg =" ")
{
    if( ( index < 0 ) || ( index > slist->entry_count-1) ) {
       if (index < 0) {
          cerr<<"List-Class-Error: "<<msg<<" : index less than zero \n";
          abort();
       }
       else {
          cerr<<"List-Class-Error: "<<msg
              <<" :index greater than number of entries = "
              << slist->entry_count<<"\n";
          abort();
       }
    }
    return( index );
}


// This defines the list base basic find operation for single link lists
// It moves through the list and returns a pointer to the element with the
// requested index.
// This routines assumes a non-empty list.

inline slink* slist_base::base_bfind(slist_base * const slist,int index)
{
    slink* ret_val = slist->last->next;
    while( index-- > 0)
       ret_val = ret_val->next;
    return ret_val;
}





// This defines the list base find operation, which returns a pointer to the
// element with the requested index (default = 0) or returns NULL, if there
// is none.

inline slink* slist_base::base_find(slist_base * const slist, int index = 0)
{
    if( slist->last == 0 )
       return 0;
    base_index_check(slist,index,"find");
    return base_bfind(slist,index);
}


// This defines the list base get operations, which returns a pointer to the
// element with the requested index (default = 0) and removes it from the
// list. If there is no element, NULL is returned.
//
// Note: Do not forget to delete the elements explicitely, after removing
//       them by the get operation !!!!!!!

inline slink* slist_base::base_get(slist_base * const slist, int index = 0,
                                                             char* msg="get")
{
    if( slist->last == 0 )
       return 0;
    base_index_check(slist,index,msg);
    if( index == 0 ) {
       slink* ret_val  = base_bfind(slist, 0 );
       if ( ret_val == slist->last )  // remove last and only element
          slist->last = 0;
       else                        // remove first element
          slist->last->next = ret_val->next;
       slist->entry_count--;
       return ret_val;
    }
    slink* ret_prev = base_bfind(slist,index - 1);
    slink* ret_val  = ret_prev->next;
    ret_prev->next = ret_val->next;
    if( ret_val == slist->last )
       slist->last=ret_prev;
    slist->entry_count--;
    return ret_val;
}


// This defines the list base clear operation. It clears the list elements
// from the list, but the objects will still exist after this routine is
// finished. The list element objects are not cleared.

inline void slist_base::base_clear(slist_base * const slist) 
{
  slist->last=0;
    slist->entry_count=0;
}


// This defines the list base next operation. If the list is empty or the
// current position is at the end of the list, NULL is returned. If an
// optional flag is set the current position is reset to the beginning.

inline slink* slist_base::base_next(slist_base * const slist, slink* curr,
                                                      ListBool reset = False)
{
    if( slist->last == 0)
        return 0;
    if( reset )
        curr = slist->last;
    else if( curr == slist->last )
        return 0;
    return curr->next;
}

// This defines the list base contains operation. If there are any elements
// matching the supplied data, True is returned. An empty list always returns
// a False.
//
// Note that this is a list item pointer comparison. Non-Intrusive lists
// must check the data stored in with the link.

inline ListBool slist_base::base_contains(slist_base * const slist, slink* datum)
{
    slink* rover = 0;
    rover = base_next(slist, rover, True);
    while( rover != 0 ) {
       if( rover == datum )
          return True;
       rover = base_next(slist, rover, False);
    }
    return False;
}

#endif

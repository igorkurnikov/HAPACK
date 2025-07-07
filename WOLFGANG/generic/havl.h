#ifndef __HAVL_HPP
#define __HAVL_HPP

//////////////////////////////////////////////////////////////////////
//
// AVL (Adelson-Velski and Landis) trees are the 'standard' computer
// science version of balanced binary trees.  This interface describes
// general purpose AVL trees
//

#include <assert.h>
#include "io_incl.h"

class Node;

template<class Id,class Node> class HAVLTree;

template<class Id,class Node> class HAVLNode
{
public:
    friend class HAVLTree<Id,Node>;

    Node* nd;

    HAVLNode<Id,Node> *&ptr(int a) { return (a < 0)? left: right; }
    HAVLNode<Id,Node> *left, *right;
    signed char balance;

    HAVLNode(Node* node) : nd(node) 
    { left = right = 0; balance = 0; }
   ~HAVLNode() { delete nd; }
};

template <class Id,class Node> class HAVLpath 
{
     friend class HAVLTree<Id,Node>;
     int level;
     HAVLNode<Id,Node> *ptr[32];
     int descent[32];

     HAVLpath() { level = 0; }
void push(HAVLNode<Id,Node>* p, int d)  
     { ptr[level] = p; descent[level] = d; ++level; }
void pop (HAVLNode<Id,Node>* &p, int &d)
     { --level; p = ptr[level]; d = descent[level]; } 
};

enum order { inorder, revorder, preorder, postorder };
enum direction { higher = +1, lower = -1 };

// really dirty, but i have no clue how else

typedef void (*vf_t)(void*, void *);

template<class Id,class Node> class HAVLTree 
{
  unsigned int no;
  HAVLNode<Id,Node> *root;
public:
      HAVLTree() { root = 0; no = 0; }
     ~HAVLTree() { cleanup(); }
int   insert(Node* newnode);
int   remove(const Id&   id);
int   size() const { return no; }

Node* lookup(const Id& id) const
{
  HAVLNode<Id,Node> *p = root;

  while(p) 
  {
    int cmp;
    if ((cmp = id.compare(*(p->nd))) == 0)
        return p -> nd;
    else
        p = p -> ptr(cmp);
  }
  return 0;
}


Node* max() const;
Node* min() const;

int   concatenate(HAVLTree<Id,Node> &);
void  traverse(void (*vf)(Node*,void*), void *supplementary, order o = inorder)
      { _trav(vf, supplementary, o, root); }
void   read (FILE* f);
void   write(FILE* f);
void   cleanup(HAVLNode<Id,Node> * = 0);

private:
       HAVLNode<Id,Node> *rotate(HAVLNode<Id,Node> *, int);
void  _trav(void (*F)(Node*,void*), void *, order o, HAVLNode<Id,Node> * = 0);
int    concat (HAVLTree<Id,Node> &a, HAVLTree<Id,Node> &b,HAVLNode<Id,Node>*);
int    depth  (HAVLNode<Id,Node> *) const;
};

/*
template <class Id,class Node> void printnode__generic(Node* node,ostream* os)
{
  (*os) << *node << endl;
}
*/
/*
template <class Id,class Node> ostream& operator<<(ostream& os,HAVLTree<Id,Node> tree)
{
  typedef void (*TRAV)(Node*,void*);
  tree.traverse(  printnode__generic,(void*) &os);
  return os;
}
*/

/*
template <class Id,class Node> void writenode__generic(Node* node,FILE* file)
{
  int cnt = fwrite((char*) node,sizeof(Node),1,file);
  cout << "write: " << node << endl;
  assert(cnt == 1);
}
*/
  /* what happens here: 
      we define 2 types, one is the pointer to a function
      which fits as the first argument for traverse.
      We cannot directly cast writenode__generic into 
      this argument, since the function has not been
      instantiated. Therefore we define a second type, 
      the pointer of a function, which has exactly
      the desired type of writenode generic. The following
      assginment intantiates the function template.
      The pointerr thus obtained is then cast onto the 
      proper type for the traverse function.
      
  */

/*
template <class Id,class Node> void HAVLTree<Id,Node>::write(FILE* file)
{
  typedef void (* F)(Node*,void*);
  typedef void (*FN)(Node*,FILE*);
  FN func = writenode__generic;
  traverse( (F) func,(void*) file);
  Node node;
  writenode__generic(&node,file);
}
*/
/*
template <class Id,class Node> void HAVLTree<Id,Node>::read(FILE* file)
{
   Node node;
   while(!feof(file))
   {
     node.read(file);

     if (!node.allowed()) break;
     insert(node);
   }
}
*/
template <class Id,class Node> void HAVLTree<Id,Node>::cleanup(HAVLNode<Id,Node> *p)
{
  if (!p) 
  {
    if (!root) return;
    p = root;
  }
  if (p->left)  cleanup(p->left);
  if (p->right) cleanup(p->right);  
  if (p == root)
    root = 0;
  delete p;
  root = 0;
  
}

#include <assert.h>

template <class Id,class Node> int HAVLTree<Id,Node>::insert(Node* node)
{
    HAVLNode<Id,Node> *t = 0, *s, *p, *q, *r;

    s = p = root;
    int a, gpd;

    // s marks the spot where rebalancing may be required
    // t is the father of s
    
    if (!root) // deal with the empty tree
    {
      root = s = q = new HAVLNode<Id,Node>(node);
      assert(root != 0);
      no++;
      return 1;
    }
      
    // find the spot to insert the new node
    for (;;) 
    {
      // compare new node with current node
      a = node -> compare(*(p -> nd));
      if (a==0) return 0; // duplicate elements, we are done

      if (a < 0) // element is smaller go to the left
      {
        q = p->left; 
        if (!q)       // no left element, make one
        {
          q = new HAVLNode<Id,Node>(node);
          assert(q != 0);
          p->left = q;
          break;
        }
      } 
      else // element is larger, go to the right
      {
        q = p->right;
        if (!q)  // no right element, make one
        { 
           q = new HAVLNode<Id,Node>(node);
           assert(q != 0);
           p->right = q;
           break;
        }
      }
      
      /* remember the last node, which has a nonzero balance factor
                 t points to the this node,
                 s points to its son in the proper direction 
         gdp is the direction 
      */
      if (q->balance) 
      {
          t = p;     // on exit t is the old node
          gpd = a;
          s = q;    //
      }
      p = q; // continue with the next node in the proper direction
    }
        
    /* 
       now deal with the balancing: 
    */
    
    // is the new node less than his father ?
    
    if(node -> compare(*(s -> nd)) < 0)  // we are on a left branch
    {
      r = p = s->left; a = -1;
    }
    else
    {
      r = p = s->right;
      a = +1;
    }
    while(p != q) 
    {
       if (node -> compare(*(p -> nd)) < 0) 
       {
          p->balance = -1;
          p = p->left;
       } 
       else 
       {
          p->balance = +1;
          p = p->right;
       }
    }
    
    if (!s->balance) 
    {
        s->balance = a;
	no++;
	return 1;
    }
    if (s->balance == -a) 
    {
        s->balance = 0;
	no++;
	return 1;
    }
    
    if (t)
        t->ptr(gpd) = rotate(s,a);
    else
        root = rotate(s,a);
    no++;
    return 1;
}

template <class Id,class Node> HAVLNode<Id,Node>* HAVLTree<Id,Node>::rotate(HAVLNode<Id,Node> *s, int sense)
{
    HAVLNode<Id,Node> *rv, *r = s->ptr(sense);
    if (r->balance == 0) { // single rotation, ends up unbalanced
        rv = r;
        s->ptr(sense) = r->ptr(-sense);
        r->ptr(-sense) = s;
        s->balance = sense;
        r->balance = -sense;
    } else if (r->balance == sense) { // single rotation
        rv = r;
        s->ptr(sense) = r->ptr(-sense);
        r->ptr(-sense) = s;
        s->balance = r->balance = 0;
    } else {                // double rotation
        rv = r->ptr(-sense);
        r->ptr(-sense) = rv->ptr(sense);
        rv->ptr(sense) = r;
        s->ptr(sense) = rv->ptr(-sense);
        rv->ptr(-sense) = s;
        if (!rv->balance) {
            s->balance = 0;
            r->balance = 0;
        } else if (rv->balance == sense) {
            s->balance = -sense;
            r->balance = 0;
        } else {
            s->balance = 0;
            r->balance = sense;
        }
        rv->balance = 0;
    }
    return rv;
}
/*

template <class Id,class Node> void HAVLTree<Id,Node>::dump(HAVLNode<Id,Node> *p) const
{
   if (!p) 
   {
     if (root) 
     {
        cout << "root " << (void *) root << endl;
        p = root;
     } 
     else
     {
        cout << " Tree Empty. " << endl;
        return;
     }
  }
  cout << (void *) p << ' ' <<  (int) p->balance << ' ' 
       << (void *) p->left << ' ' 
       << (void *) p->right << ' ' << (Node&) *p << endl;
  
  if (p->left) dump(p->left);
  if (p->right) dump(p->right);
}

*/
/*
template <class Id,class Node> Node* HAVLTree<Id,Node>::lookup(Node& v) const
{
  HAVLNode<Id,Node> *p = root;

  while(p) 
  {
    int cmp;
    if ((cmp = v.compare(*p)) == 0)
        return p;
    else
        p = p->ptr(cmp);
  }
  return 0;
}
*/
template<class Id,class Node> int HAVLTree<Id,Node>::remove(const Id& v)
{    
    HAVLNode<Id,Node> *p = 0, *q = root, *r, *s, *t;

    int pd, a, tlev;
    HAVLpath<Id,Node> path;
    for (;;) {
        if (q == 0)
            return 0;       // no such item
        a = v.compare(*(q -> nd));
        if (!a) break;      // found
        a = (a < 0)? -1: +1;
        path.push(q,a);     // note the path
        q = q->ptr(a);
    }
    if (path.level) {
        p = path.ptr[path.level-1];
        a = path.descent[path.level-1];
    }
    pd = a;

    if (!q->right) {
        if (p) {
            path.level--;       // already got this path entry
            p->ptr(pd) = q->left;
            r = p;
        } else {
            root = q->left;
            r = root;
            if (!r) 
	    {
                delete q;
		no--;
		return 1;
            }
        }
    } else {
        r = q->right;
        a = +1;
        if (!r->left) {
            r->left = q->left;
            r->balance = q->balance;
            if (p) {
                p->ptr(pd) = r;
            } else {
                root = r;
            }
        } else {
            path.push(q,1);         // add target node to stack
            tlev = path.level-1;
            for (s = r; s->left;) {
                path.push(s,-1);
                a = -1;
                s = s->left;
            }
            path.pop(r,a);
            s->left = q->left;
            r->left = s->right;
            s->right = q->right;
            s->balance = q->balance;
            if (p)
                p->ptr(pd) = s;
            else 
                root = s;
            path.ptr[tlev] = s;     // adjust path to include switched node
        }
    }
    for (;;) {
        int done = 0;
        t = r;
// a side subtree of r has decreased in height
        if (r->balance == a) {
            r->balance = 0;
        } else if (r->balance == 0) {
            r->balance = -a;
            delete q;
	    no--;
	    return 1;
        } else {
            done = (r->ptr(-a)->balance == 0);
            t = rotate(r,-a);
        }
        if (path.level) {
            path.pop(r,a);
            r->ptr(a) = t;
        } else {
            root = t;
            break;
        }
        if (done) break;
    }
    delete q;
    no--;
    return 1;
}

// HAVL trees are balanced, so recursive functions should not be problematic

template <class Id,class Node> void 
HAVLTree<Id,Node>::_trav(void (*function)(Node*, void*), 
			 void *supplementary,order o,HAVLNode<Id,Node> *p)
{
    if (!p) {
        if (!root)
            return;
        p = root;
    }
    
    switch (o) {
    case inorder:
        if (p->left) _trav(function, supplementary, o, p->left);
        function(p -> nd, supplementary);
        if (p->right) _trav(function, supplementary, o, p->right);
        return;
    case revorder:
        if (p->right) _trav(function, supplementary, o, p->right);
        function(p -> nd, supplementary);
        if (p->left) _trav(function, supplementary, o, p->left);
        return;
    case preorder:
        function(p -> nd, supplementary);
        if (p->left) _trav(function, supplementary, o, p->left);
        if (p->right) _trav(function, supplementary, o, p->right);
        return;
    case postorder:
        if (p->right) _trav(function, supplementary, o, p->right);
        function(p -> nd, supplementary);
        if (p->left) _trav(function, supplementary, o, p->left);
        return;
    }
}
   

#endif

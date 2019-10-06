#include <set>
#include <cassert>
#include <iterator>
#include <iostream>
#include <common.hpp>
#include <objfun.hpp>
#include "morf.hpp"
#include "atmv.hpp"
#include "iter.hpp"
#include "rang.hpp"

namespace mort
{
    using namespace boost;

    atmvec::atmvec()
    {
    }

    atmvec::atmvec(int size, const atom_t& a0) :
        vector<atom_t> (size, a0)
    {
        assert( a0.cmpid() == ATOM );
    }

    atmvec::atmvec(const atom_range& range) :
        vector<atom_t> (range.begin(), range.end())
    {
    }

    atmvec::atmvec(const morf_t& o)
    {
        if (o.cmpid() == ATOM)
        {
            push_back(o);
        }
        else
        {
            push_atom(o);
        }
    }

    atmvec::atmvec(const atom_t& a0, const atom_t& a1) :
        vector<atom_t> (2, a0)
    {
        at(1) = a1;
    }

    atmvec::atmvec(const atom_t& a0, const atom_t& a1,
            const atom_t& a2) :
        vector<atom_t> (3, a0)
    {
        at(1) = a1;
        at(2) = a2;
    }

    atmvec::atmvec(const atom_t& a0, const atom_t& a1,
            const atom_t& a2, const atom_t& a3) :
        vector<atom_t> (4, a0)
    {
        at(1) = a1;
        at(2) = a2;
        at(3) = a3;
    }

    atmvec::atmvec(const atom_t& a0, const atom_t& a1,
            const atom_t& a2, const atom_t& a3, const atom_t& a4) :
        vector<atom_t> (5, a0)
    {
        at(1) = a1;
        at(2) = a2;
        at(3) = a3;
        at(4) = a4;
    }

    atmvec::atmvec(const atom_t& a0, const atom_t& a1,
            const atom_t& a2, const atom_t& a3, const atom_t& a4,
            const atom_t& a5) :
        vector<atom_t> (6, a0)
    {
        at(1) = a1;
        at(2) = a2;
        at(3) = a3;
        at(4) = a4;
        at(5) = a5;
    }

    atmvec::~atmvec()
    {
    }

    atmvec& atmvec::operator=(const atmvec& rhs)
    {
        if (&rhs != this)
        {
            atmvec tmp(rhs);
            root_t::swap(tmp);
            vector<atom_t>::swap(tmp);
        }

        return *this;
    }

    int atmvec::count(const atom_t& atom) const
    {
        return std::count(begin(), end(), atom);
    }

    void atmvec::push_atom(const morf_t& o)
    {
        if( o.cmpid()==ATOM )
        {
            push_back( o );
        }
        else
        {
            insert(end(), o.atom_begin(), o.atom_end());
        }
    }

    std::ostream& operator<<(std::ostream& os, const atmvec& av)
    {
        for(unsigned int i = 0; i < av.size(); ++i)
        {
            os << av[i].name() << " ";
        }

        os << "\n";
        return os;
    }

    bool atmvec::is_ring() const
    {
        return size() > 2 && front().is_connected_to(back());
    }

    bool atmvec::is_arom() const
    {
        const_iterator ai = begin();
        for (; ai != end() && is_sp2(*ai); ++ai)
            ;

        if (ai != end())
            return false;

        int total_pie = sum(begin(), end(), 0, get_atom_pie);

        return (total_pie % 4) == 2;
    }

    bool find_path( const atom_t& a, const atom_t& end, atmvec& path, std::set<int>& visited)
    {
        if (visited.count(a.absid()))
        {
            return false;
        }

        visited.insert(a.absid());

        path.push_back(a);

        if (a == end)
        {
            return true;
        }

        atomiter_t n = a.atom_begin();
        for(; n != a.atom_end(); ++n)
        {
            if(find_path(*n, end, path, visited))
            {
                return true;
            }
        }

        path.pop_back();
        return false;
    }

    atmvec find_path(const atom_t& start, const atom_t& dest)
    {
        std::set<int> visited;

        atmvec path;

        return find_path(start, dest, path, visited) ? path : atmvec();
    }

} // namespace mort



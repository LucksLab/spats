
#ifndef __SPATS_LOOKUP_HPP_INCLUDED__
#define __SPATS_LOOKUP_HPP_INCLUDED__

#include "target.hpp"


class Lookup :
    public virtual FragmentStorage
{
private:
    int m_expectedSize;
    int m_modulus;
    int m_listSize;
    int m_count;
    FragmentResult ** m_table;
    void rebuild(int minSize);

public:
    Lookup(int expectedSize) :
        m_expectedSize(expectedSize), m_table(NULL), m_count(0), m_listSize(0)
    {
        rebuild(0);
    }

    virtual ~Lookup()
    {
        delete [] m_table;
    }

    int count() const { return m_count; }

    bool insert(FragmentResult * f);
    FragmentResult * find(Fragment *f, int start_index = 0) const;
    void dump() const;
};


#endif // __SPATS_R1TREE_HPP_INCLUDED__


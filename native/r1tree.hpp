
#ifndef __SPATS_R1TREE_HPP_INCLUDED__
#define __SPATS_R1TREE_HPP_INCLUDED__

#include "target.hpp"

class R1Tree;


class R1Tree :
    public virtual FragmentStorage
{
private:

    Fragment m_prefix;
    int m_r1Length;
    int m_prefixLength;
    R1Tree * m_next[4];
    FragmentResult * m_leaf;

    R1Tree(R1Tree * parent, uint64_t nextNt, FragmentResult * f);
    R1Tree(R1Tree * parent, FragmentResult * f1, FragmentResult * f2, int match_length);

    inline void
    clearNext()
    {
        for (int i = 0; i < 4; ++i)
            m_next[i] = NULL;
    }

public:
    R1Tree(int r1Length) :
        m_r1Length(r1Length), m_prefixLength(0), m_leaf(NULL)
    {
        clearNext();
    }

    virtual ~R1Tree()
    {
        for (int i = 0; i < 4; ++i)
            if (m_next[i])
                delete m_next[i];
    }
    int count() const { return 0; }
    bool insert(FragmentResult * f);
    FragmentResult * find(Fragment * f, int start_index = 0) const;
    void dump() const;
};


#endif // __SPATS_R1TREE_HPP_INCLUDED__


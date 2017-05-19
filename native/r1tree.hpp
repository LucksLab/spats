
#ifndef __SPATS_R1TREE_HPP_INCLUDED__
#define __SPATS_R1TREE_HPP_INCLUDED__

#include "seq.hpp"

class Target
{
public:
    int m_identifier;
    std::string m_name;
    Target(int id, std::string name) : m_identifier(id), m_name(name) {}
};

class FragmentResult
{
public:
    Fragment m_fragment;
    Target * m_target;
    int m_site;
};

class R1Tree;

class R1Tree
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

    ~R1Tree()
    {
        for (int i = 0; i < 4; ++i)
            if (m_next[i])
                delete m_next[i];
    }

    bool insert(FragmentResult * f);
    FragmentResult * find(Fragment &f) const;
    void dump() const;
};


#endif // __SPATS_R1TREE_HPP_INCLUDED__


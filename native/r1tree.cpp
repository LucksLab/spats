
#include <stdio.h>

#include "ats.hpp"
#include "r1tree.hpp"


R1Tree::R1Tree(R1Tree * parent, uint64_t nextNt, FragmentResult * f) :
    m_r1Length(parent->m_r1Length), m_prefixLength(1 + parent->m_prefixLength), m_leaf(f)
{
    clearNext();
    m_prefix.clone(&f->m_fragment);
}

R1Tree::R1Tree(R1Tree * parent, FragmentResult * f1, FragmentResult * f2, int match_length) :
    m_r1Length(parent->m_r1Length), m_prefixLength(match_length), m_leaf(NULL)
{
    clearNext();
    m_prefix.clone(&f1->m_fragment);
    uint64_t nextNt1 = f1->m_fragment.at(m_prefixLength);
    uint64_t nextNt2 = f2->m_fragment.at(m_prefixLength);
    ATS_ASSERT(nextNt2 != nextNt2);
    m_next[nextNt1] = new R1Tree(this, nextNt1, f1);
    m_next[nextNt2] = new R1Tree(this, nextNt2, f2);
}

bool
R1Tree::insert(FragmentResult * f)
{
    if (!m_prefix.equals(&f->m_fragment, m_prefixLength)) {
        int length = m_prefixLength;
        while (length > 0  &&  !m_prefix.equals(&f->m_fragment, length))
            --length;
        R1Tree * newNode = new R1Tree(m_r1Length);
        newNode->m_prefix.clone(&m_prefix);
        newNode->m_prefixLength = m_prefixLength;
        for (int i = 0; i < 4; ++i) {
            newNode->m_next[i] = m_next[i];
            m_next[i] = NULL;
        }
        m_prefixLength = length;
        m_next[m_prefix.at(m_prefixLength)] = newNode;
        newNode->m_leaf = m_leaf;
        m_leaf = NULL;
        /* fall-through will end up in last case below, just inserting new node for f */
    }
    uint64_t nextNt = f->m_fragment.at(m_prefixLength);
    if (m_next[nextNt])
        return m_next[nextNt]->insert(f);
    if (m_leaf) {
        uint64_t leafNextNt = m_leaf->m_fragment.at(m_prefixLength);
        if (leafNextNt == nextNt) {
            int length = m_prefixLength + 1;
            while (m_leaf->m_fragment.equals(&f->m_fragment, length))
                ++length;
            R1Tree * newNode = new R1Tree(this, m_leaf, f, length);
            m_next[nextNt] = newNode;
            m_leaf = NULL;
        }
        else {
            R1Tree * newLeafNode = new R1Tree(this, leafNextNt, m_leaf);
            m_next[leafNextNt] = newLeafNode;
            R1Tree * newNode = new R1Tree(this, nextNt, f);
            m_next[nextNt] = newNode;
            m_leaf = NULL;
        }
    }
    else {
        R1Tree * newNode = new R1Tree(this, nextNt, f);
        m_next[nextNt] = newNode;
    }
    return true;
}

FragmentResult *
R1Tree::find(Fragment * f, int start_index) const
{
    if (!m_prefix.equals(f, m_prefixLength, start_index))
        return NULL;
    uint64_t nextNt = f->at(m_prefixLength);
    if (m_next[nextNt])
        return m_next[nextNt]->find(f, m_prefixLength + 1);
    else if (m_leaf  &&  f->equals(&m_leaf->m_fragment, m_r1Length))
        return m_leaf;
    else
        return NULL;
}

void
R1Tree::dump() const
{
    assert(0); // TODO
#if 0
    std::string indent(m_prefixLength, ' ');
    printf("%spfx: %s\n", indent.c_str(), m_prefix.string(m_prefixLength).c_str());
    if (m_leaf)
        printf("%sLEAF: %s\n", indent.c_str(), m_leaf->m_fragment.string(m_r1Length).c_str());
    for (int i = 0; i < 4; ++i) {
        if (m_next[i]) {
            printf("%s%c:\n", indent.c_str(), nt_bits_to_ch(i));
            m_next[i]->dump();
        }
    }
#endif
}


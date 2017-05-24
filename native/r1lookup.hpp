
#ifndef __SPATS_R1LOOKUP_HPP_INCLUDED__
#define __SPATS_R1LOOKUP_HPP_INCLUDED__

#include "target.hpp"


class R1Lookup
{
private:
    Targets * m_targets;
    std::string m_adapter_b;
    int m_r1_length;
    int m_allowed_errors;
    FragmentStorage * m_storage;

    void build();

public:

    R1Lookup(Targets * targets, std::string adapter_b, int r1_length, int allowed_errors = 0):
        m_targets(targets), m_adapter_b(adapter_b), m_r1_length(r1_length), m_allowed_errors(allowed_errors), m_storage(NULL)
    {
        build();
    }

    ~R1Lookup()
    {
        if (NULL != m_storage)
            delete m_storage;
    }

    FragmentResult * find(Fragment * f) const;
    void dump(void) const;
};


class R2Lookup
{
private:
    int m_r2_length;
    int m_allowed_errors;
    int m_numTargets;
    FragmentStorage ** m_targetStores;

public:

    R2Lookup(int r2_length, int allowed_errors = 0):
        m_r2_length(r2_length), m_allowed_errors(allowed_errors), m_targetStores(NULL), m_numTargets(0)
    {
    }

    ~R2Lookup()
    {
        for (int i = 0; i < m_numTargets; ++i)
            delete m_targetStores[i];
        delete [] m_targetStores;
    }

    void addTarget(Target * t);
    FragmentResult * find(Fragment * f, Target * t) const;
    void dump(void) const;
};


#endif // __SPATS_R1TREE_HPP_INCLUDED__


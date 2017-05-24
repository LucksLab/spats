
#ifndef __SPATS_R1LOOKUP_HPP_INCLUDED__
#define __SPATS_R1LOOKUP_HPP_INCLUDED__

#include "target.hpp"


class R1Lookup
{
private:
    Target * m_target;
    std::string m_adapter_b;
    int m_r1_length;
    int m_allowed_errors;
    FragmentStorage * m_storage;

    void build();

public:

    R1Lookup(Target * target, std::string adapter_b, int r1_length, int allowed_errors = 0):
        m_target(target), m_adapter_b(adapter_b), m_r1_length(r1_length), m_allowed_errors(allowed_errors), m_storage(NULL)
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


#endif // __SPATS_R1TREE_HPP_INCLUDED__


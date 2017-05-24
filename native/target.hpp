
#ifndef __SPATS_TARGET_HPP_INCLUDED__
#define __SPATS_TARGET_HPP_INCLUDED__

#include "seq.hpp"

class Target
{
    int m_identifier;
    std::string m_name;
    std::string m_seq;

public:
    Target(int id, std::string name, std::string seq) : m_identifier(id), m_name(name), m_seq(seq) {}
    int identifier() const { return m_identifier; }
    const std::string& name() const { return m_name; }
    const std::string& seq() const { return m_seq; }
    int n() const { return m_seq.size(); }
    std::string subseq(int start, int length) const { return m_seq.substr(start, length); }
};

class Targets
{
    int m_numTargets;
    Target ** m_targets;
public:
    Targets() : m_numTargets(0), m_targets(NULL) { }
    ~Targets()
    {
        for (int i = 0; i < m_numTargets; ++i)
            delete m_targets[i];
        if (m_targets)
            delete [] m_targets;
    }
    int size() const { return m_numTargets; }
    void addTarget(Target * target);
    Target * target(int idx) const { return m_targets[idx]; }
    void parse(const char * fasta_path);
};

class FragmentResult
{
public:
    Fragment m_fragment;
    Target * m_target;
    int m_site;
    int m_errors;
    int m_count;
    FragmentResult() { }
    FragmentResult(Target * target, const char * text, int length, int site, int errors = 0) :
        m_target(target), m_fragment(text, length), m_site(site), m_errors(errors), m_count(0) { }
    FragmentResult(FragmentResult * other, int errors = -1) :
        m_target(other->m_target), m_site(other->m_site), m_errors(-1 == errors ? other->m_errors : errors), m_count(0)
    {
        m_fragment.clone(&other->m_fragment);
    }
};

class FragmentStorage
{
public:
    virtual int count() const = 0;
    virtual bool insert(FragmentResult * f) = 0;
    virtual FragmentResult * find(Fragment * f, int start_index = 0) const = 0;
    virtual void dump() const = 0;
    virtual ~FragmentStorage() { }
};

std::string
reverse_complement(std::string seq);

#endif // __SPATS_TARGET_HPP_INCLUDED__

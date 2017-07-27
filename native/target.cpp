
#include "target.hpp"
#include "parse.hpp"
#include "ats.hpp"

#define BUFFER_SIZE 2048

std::string
reverse_complement(std::string seq)
{
    char buffer[BUFFER_SIZE];
    int len = seq.size();
    assert(len < BUFFER_SIZE);
    const char * cstr = seq.c_str();
    char ch = 0;
    for (int i = 0; i < len; ++i) {
        ch = cstr[len - (i + 1)];
        switch (ch) {
        case 'A': ch = 'T'; break;
        case 'C': ch = 'G'; break;
        case 'G': ch = 'C'; break;
        case 'T': ch = 'A'; break;
        default: break;
        }
        buffer[i] = ch;
    }
    return std::string(buffer, len);
}

Targets * g_parsing_targets = NULL;

bool
Targets_fa_handler(const char * name, const char * seq)
{
    g_parsing_targets->addTarget(new Target(g_parsing_targets->size(), name, seq));
    return true;
}

void
Targets::addTarget(Target * target)
{
    Target ** targets = new Target *[1 + m_numTargets];
    for (int i = 0; i < m_numTargets; ++i)
        targets[i] = m_targets[i];
    if (m_targets)
        delete [] m_targets;
    m_targets = targets;
    targets[m_numTargets] = target;
    ++m_numTargets;
}

void
Targets::parse(const char * fasta_path)
{
    g_parsing_targets = this;
    fasta_parse(fasta_path, &Targets_fa_handler);
}

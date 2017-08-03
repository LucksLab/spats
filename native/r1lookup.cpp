
#include "r1lookup.hpp"

#include "lookup.hpp"
#include "ats.hpp"


void
R1Lookup::build()
{
    int length = m_r1_length;
    int expected_size = length + (m_allowed_errors * length * 9);
    m_storage = new Lookup(expected_size * m_targets->size());
    for (int target_idx = 0; target_idx < m_targets->size(); ++target_idx) {
        Target * target = m_targets->target(target_idx);
        int tlen = target->n();
        Fragment f;
        std::string rc = reverse_complement(target->seq());
        for (int i = 0; i < length + 1; ++i) {
            std::string candidate = rc.substr(0, i) + m_adapter_b.substr(0, length - i);
            ATS_VERBOSE(" C %d: %s", i, candidate.c_str());
            FragmentResult * fr = new FragmentResult(target, candidate.c_str(), length, i == length ? -1 : tlen - i);
            m_storage->insert(fr);
            if (m_allowed_errors > 0) {
                assert(1 == m_allowed_errors);
                for (int j = 0; j < m_r1_length; ++j) {
                    uint64_t cur = fr->m_fragment.at(j);
                    for (int k = 0; k < 4; ++k) {
                        FragmentResult * ins = new FragmentResult(fr, 1);
                        ins->m_fragment.insert(j, k);
                        m_storage->insert(ins);
                        if (cur == k)
                            continue;
                        FragmentResult * toggle = new FragmentResult(fr, 1);
                        toggle->m_fragment.set(j, k);
                        m_storage->insert(toggle);
                    }
                    FragmentResult * del = new FragmentResult(fr, 1);
                    del->m_fragment.del(j);
                    m_storage->insert(del);
                }
            }
        }
    }
    ATS_DEBUG("Lookup table: %d R1 entries...\n", m_storage->count());
    //m_storage->dump();
}

void
R1Lookup::build_cotrans()
{
    int minimum_target_length = 3;
    int linker_len = m_linker->len();
    Target * target = m_targets->target(0);
    int tlen = target->n();
    int r1_match_len = m_r1_length - 4;
    int target_left = r1_match_len - linker_len - minimum_target_length + 1;
    int expected_size = tlen * target_left;
    m_storage = new Lookup(expected_size, true);

    for (int end = minimum_target_length; end < tlen + 1; ++end) {
        std::string target_subseq(target->seq().substr(0, end));
        Fragment f;
        for (int i = 0; i < target_left; ++i) {
            int tstart = i - (r1_match_len - linker_len);
            if (tstart + end < 0)
                continue;
            std::string r1_rc_match(target_subseq.substr(end + tstart, 0 - tstart) + m_linker->string());
            std::string r1_match(reverse_complement(r1_rc_match) + m_adapter_b.substr(0, i));
            ATS_VERBOSE(" C %d/%d: %s", end, i, r1_match.c_str());
            FragmentResult * fr = new FragmentResult(target, r1_match.c_str(), r1_match_len, end, 0, i);
            m_storage->insert(fr);
        }
    }

    ATS_DEBUG("Lookup table: %d R1 entries...\n", m_storage->count());
    //m_storage->dump();
}

FragmentResult * 
R1Lookup::find(Fragment *f) const
{
    return m_storage->find(f);
}

void
R1Lookup::dump(void) const
{
    m_storage->dump();
}


void
R2Lookup::addTarget(Target * t)
{
    assert(t->identifier() == m_numTargets);
    FragmentStorage ** stores = new FragmentStorage *[1 + m_numTargets];
    for (int i = 0; i < m_numTargets; ++i)
        stores[i] = m_targetStores[i];
    if (m_targetStores)
        delete [] m_targetStores;
    m_targetStores = stores;

    int n = t->n();
    std::string tseq = t->seq();
    int expected_size = n + (m_allowed_errors * n * 9);
    FragmentStorage * store = new Lookup(expected_size);
    std::string candidate;

    for (int i = 0; i < n; ++i) {
        if (i + m_r2_length <= n)
            candidate = tseq.substr(i, m_r2_length);
        else
            candidate = tseq.substr(i);
        FragmentResult * fr = new FragmentResult(t, candidate.c_str(), -1, i);
        store->insert(fr);
    }

    m_targetStores[m_numTargets] = store;
    ++m_numTargets;
}

FragmentResult * 
R2Lookup::find(Fragment * f, Target * t) const
{
    return m_targetStores[t->identifier()]->find(f);
}

void
R2Lookup::dump(void) const
{
}

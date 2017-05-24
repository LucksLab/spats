
#include "r1lookup.hpp"

#include "r1tree.hpp"
#include "lookup.hpp"
#include "ats.hpp"


void
R1Lookup::build()
{
    int expected_size = m_r1_length + (m_allowed_errors * m_r1_length * 9);
    m_storage = new Lookup(expected_size);
    int length = m_r1_length;
    int tlen = m_target->n();
    Fragment f;
    ATS_VERBOSE(" L: %d", length);
    std::string rc = reverse_complement(m_target->seq());
    for (int i = 0; i < length + 1; ++i) {
        std::string candidate = rc.substr(0, i) + m_adapter_b.substr(0, length - i);
        ATS_VERBOSE(" C %d: %s", i, candidate.c_str());
        FragmentResult * fr = new FragmentResult(m_target, candidate.c_str(), length, i == length ? -1 : tlen - i);
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
    printf("Lookup table: %d R1 entries...\n", m_storage->count());
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

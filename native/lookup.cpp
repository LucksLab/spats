
#include "lookup.hpp"
#include "ats.hpp"
#include <stdlib.h>


void
Lookup::rebuild(int minSize)
{
    FragmentResult ** originalTable = m_table;
    int originalSize = m_listSize;

    int maxSize = std::max(std::max(m_count, m_expectedSize), 4);
    maxSize = maxSize << 3; // give a little more room

    m_modulus = 1;
    while (m_modulus < maxSize)
        m_modulus = m_modulus << 1;
    m_modulus -= 1;

    m_listSize = std::max(std::max(minSize, m_listSize), 16);

    size_t size = m_modulus * m_listSize * sizeof(FragmentResult *);
    m_table = (FragmentResult **)malloc(size);
    memset(m_table, 0, size);

    if (NULL == originalTable) {
        ATS_DEBUG("table size: %dK", (int)(size >> 10));
        return;
    }

    FragmentResult * old = NULL;
    int n = 0;
    for (int m = 0; m < m_modulus; ++m) {
        int oldBase = m * originalSize;
        int curBase = m * m_listSize;
        for (int i = 0; i < originalSize; ++i) {
            old = originalTable[oldBase + i];
            if (!old)
                break;
            m_table[curBase + i] = old;
            ++n;
        }
    }
    ATS_DEBUG("********** rebuild: %d entries (new table size %dK)", n, (int)(size >> 10));

    free(originalTable);
}

bool
Lookup::insert(FragmentResult * f)
{
    uint64_t key = f->m_fragment.key();
    uint64_t rem = key % m_modulus;
    uint64_t base = rem * m_listSize;
    uint64_t idx = 0;
    FragmentResult ** entry = &m_table[base];
    while (NULL != *entry  &&  idx < m_listSize) {
        if ((*entry)->m_fragment.equals(&f->m_fragment)) {
            (*entry)->m_site = -1;
            ATS_DEBUG("dropping dup...");
            return false;
        }
        ++idx;
        entry = &m_table[base + idx];
    }
    if (idx == m_listSize) {
        rebuild(m_listSize << 1);
        return this->insert(f);
    }
    *entry = f;
    ++m_count;
    ATS_VERBOSE(" I: %s [0x%x]", f->m_fragment.string().c_str(), rem);
    return true;
}

FragmentResult * 
Lookup::find(Fragment *f, int start_index) const
{
    uint64_t key = f->key();
    uint64_t rem = key % m_modulus;
    ATS_VERBOSE("F: 0x%x", (int)rem);
    uint64_t base = rem * m_listSize;
    uint64_t idx = 0;
    FragmentResult * entry = m_table[base];
    while (NULL != entry) {
        if (f->equals(&entry->m_fragment))
            return entry;
        ++idx;
        entry = m_table[base + idx];
    }
    return NULL;
}

void
Lookup::dump() const
{
    FragmentResult * entry = NULL;
    int entryCount = 0;
    for (int i = 0; i < m_modulus; ++i) {
        entry = m_table[i * m_listSize];
        entryCount = 0;
        while (entry) {
            //printf("  %s : %d\n", entry->m_fragment.string().c_str(), entry->m_count);
            ++entryCount;
            entry = m_table[i * m_listSize + entryCount];
        }
        if (entryCount > 0)
            printf(" 0x%x : %d\n", i, entryCount);
    }
}

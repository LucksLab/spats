
#include <stdio.h>

#include "spats.hpp"
#include "parse.hpp"
#include "mask.hpp"



bool
Spats::spats_handler(Fragment * r1, Fragment * r2, const char * handle, Counters * counters,  Case * caseinfo)
{
    //printf("R1: %s\nR2: %s\nH: %s\n", r1->string().c_str(), r2->string().c_str(), handle);

    INCR_COUNTER(counters,TOTAL);

    //if (0 == COUNTER(TOTAL) % 1000000) {
    //    printf(".");
    //    fflush(stdout);
    //}
    //if (g_total > 500) {
    //    printf("%s.%s / %s...\n", handle, r1->string().c_str(), r2->string().c_str());
    //}
    //    exit(-1);

    FragmentResult * res = m_r1l->find(r1);
    if (!res) {
        //pair.failure = Failures.nomatch
        return true;
    }
    if (res->m_site < m_cotrans_minimum_length) {
        //pair.failure = Failures.cotrans_min
        return true;
    }

    if (NULL != res->m_next) {
        FragmentResult * cur = res;
        while (cur != NULL) {
            FragmentResult * other = res;
            while (other != NULL) {
                if (other != cur  &&  other->m_trim == cur->m_trim) {
                    //pair.failure = Failures.multiple_R1
                    return true;
                }
                other = other->m_next;
            }
            cur = cur->m_next;
        }
    }

    while (NULL != res) {
        if (try_lookup_hit(res, r1, r2, handle, counters, caseinfo)) {
            INCR_COUNTER(counters,MATCHED);
            return true;
        }
        res = res->m_next;
    }

    return true;
}

bool
Spats::try_lookup_hit(FragmentResult * res, Fragment * r1, Fragment * r2, const char * handle, Counters * counters, Case * caseinfo)
{
    int linker_len = m_linker->len();
    int pair_len = r2->len();
    int L = res->m_site;
    int trim = res->m_trim;
    int site = -1;
    std::string r2_seq(r2->string());
    std::string tseq(m_cotrans_target->seq());

    if (0 == trim) {
        int r2_match_len = m_r2l->r2_length();
        Fragment r2f(r2_seq.c_str(), r2_match_len);
        FragmentResult * r2_res = m_r2l->find(&r2f, m_cotrans_target);
        if (NULL != r2_res) {
            site = r2_res->m_site;
        }
        else {
            //pair.failure = Failures.nomatch
            return false;
        }
    }
    else {
        site = L - (pair_len - linker_len - 4) + trim;
    }

    // now need to verify R2 contents
    // R2: [target][linker][handle][adapter]

    int target_match_len = std::min(pair_len, L - site);
    if (target_match_len <= 0  ||  r2_seq.substr(0, target_match_len) != tseq.substr(site, target_match_len)) {
        //pair.failure = Failures.match_errors
        return false;
    }

    if (target_match_len < pair_len) {
        int linker_match_len = std::min(linker_len, pair_len - target_match_len);
        if (r2_seq.substr(target_match_len, linker_match_len) != m_linker->string().substr(0, linker_match_len)) {
            //pair.failure = Failures.linker
            return false;
        }

        if (trim > 0  &&  r2_seq.substr(r2_seq.length() - trim, trim) != m_adapter_t_rc.substr(0, trim)) {
            //pair.failure = Failures.adapter_trim
            return false;
        }
        if (pair_len - target_match_len - linker_len - 4 > trim) {
            //pair.failure = Failures.adapter_trim
            return false;
        }
    }

    if (r1->has_errors() || r2->has_errors()) {
        INCR_COUNTER(counters,INDETERMINATE);
        return false;
    }

    int mask = match_mask(handle);
    if (MASK_NO_MATCH == mask) {
        INCR_COUNTER(counters,MASK_FAILURE);
        return false;
    }

    counters->register_site(mask, L, site);
    if (NULL != caseinfo) {
        caseinfo->mask = mask;
        caseinfo->L = L;
        caseinfo->site = site;
    }

    return true;
}

void 
Spats::setup()
{
    if (NULL != m_r1l)
        return;

    ATS_ASSERT(this->m_cotrans);
    ATS_ASSERT(NULL == m_r1l  &&  NULL == m_counters);
    ATS_ASSERT(NULL != m_cotrans_target);

    m_r1l = new R1Lookup(&m_targets, m_adapter_b.c_str(), m_pair_len, 0, m_linker);
    m_r2l = new R2Lookup(12, 0);
    m_r2l->addTarget(m_cotrans_target);
    m_counters = new Counters(m_cotrans_target->n());
}

void 
Spats::run_fastq(const char * r1_path, const char * r2_path)
{
    ATS_ASSERT(NULL == m_r1l  &&  NULL == m_counters);
    setup();
    fastq_parse_spats(r1_path, r2_path, this);
}

void 
Spats::run_db(const char * db_path)
{
    ATS_ASSERT_NOT_REACHED();
}

void
Spats::run_case(Case * c)
{
    setup();
    c->mask = MASK_NO_MATCH;
    c->L = c->site = -1;
    this->spats_handler(&c->r1, &c->r2, c->handle.c_str(), m_counters, c);
}

void
Spats::addTargets(const char * path)
{
    m_targets.parse(path);
    if (m_cotrans)
        m_cotrans_target = m_targets.target(0);
}

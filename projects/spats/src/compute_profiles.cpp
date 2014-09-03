/*
 *  compute_profiles.cpp
 *  spats
 *
 *  Created by Cole Trapnell on 4/15/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#endif

#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <functional>
#include <numeric>
#include <cmath>

#include "common.h"
#include "reads.h"
#include "hits.h"
#include "bundles.h"


void print_usage()
{
	//NOTE: SPACES ONLY, bozo
	fprintf(stderr, "compute_profiles v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "Usage:   compute_profiles [options] <targets.fa> <treated_hits.sam> <untreated_hits.sam>\n");
}

struct Adducts
{
    Adducts() {}
    
    Adducts(size_t seq_len) 
    {
        adduct_counts = vector<int>(seq_len, 0);
        reactivity = vector<double>(seq_len, 0);
    }
    
    const vector<int> adducts() const { return adduct_counts; }
    
    int total_fragments() const
    {
        return accumulate(adduct_counts.begin(), adduct_counts.end(), 0);
    }
    
    int total_adducts() const
    {
        if (adduct_counts.empty())
            return 0;
        
        return accumulate((adduct_counts.begin()+1), adduct_counts.end(), 0);
    }
    
    bool register_fragment(const MateHit& fragment)
    {
        int adduct_site = fragment.left(); // fragment is 0-index, adduct array is 1-indexed, so we don't actually need to subtract
        
        //JBL - example of why don't need to subtract
        // SEQUENCE GGACU ...
        // FRAGMENT 01234 ...
        // ADDUCT   12345 ...
        // OBSERVE   GACU ...
        // fragment.left() = 1
        // increment adduct[1]++, which actually means FRAGMENT[0] which means the 5' end
        // since fragment is 0-indexed and adduct array is 1-indexed, the shift to the left is taken care of

        if (adduct_site >= 0 && adduct_site < (int)adduct_counts.size())
        {
            //JBL - only register this fragment if the left read is within the sequence and the right read aligns with the end of
            //the RNA (or RNA subsequence)
            //fragment.right() should be N-1 for an RNA subsequence of length N because of the 0-indexing
            //adduct_counts.size() should be N
            int end = min((int)adduct_counts.size(), fragment.right());
            if (end == (int)adduct_counts.size() - 1)
            {
                adduct_counts[adduct_site]++;
                return true;
            }
        }
        return false;
    }
    
private:
    vector<int> adduct_counts;
    vector<double> reactivity;
};

//code to compute reactivities starting from specific priming site
struct TargetProfile
{
    TargetProfile() {}
    
    TargetProfile(const string& name, const string& seq) : 
        _name(name), 
        _seq(seq),
        _profile_len(_seq.length() + 1),
        _treated(_profile_len),
        _untreated(_profile_len),
        _thetas(_seq.length() + 1),
        //_normalized_thetas(_seq.length()),
        _betas(_seq.length() + 1),
        _c(0.0) {}
    
    const string& name() const { return _name; }
    const string& seq() const { return _seq; }
    
    const Adducts treated() const { return _treated; }
    const Adducts untreated() const { return _untreated; }
    
    const vector<double>& thetas() const { return _thetas; }
    const vector<double>& betas() const { return _betas; }
    
    const double c() const { return _c; }
    
    bool register_fragment(const MateHit& fragment, bool treated)
    {
        if (treated)
        {
            return _treated.register_fragment(fragment);
        }
        else
        {
            return _untreated.register_fragment(fragment);
        }
    }
    
    void calc_poisson_reactivities()
    {
        
        // TargetProfile tracks an RNA of length n. arrays have n+1 entries, 
        // with index 1 corresponding to the 5'-most base, and index n 
        // corresponding to the 3'-most base.  Index 0 stores information about 
        // unmodified RNAs, where RT has fallen off the 5' end of the 
        // strand.  This convention differs from the paper, where we refer to 
        // the 5'-most base as index n, and the 3'-most base as index 1.
        if (_profile_len == 0)
            return;
        
        const vector<int>& treated_adducts = _treated.adducts();
        const vector<int>& untreated_adducts = _untreated.adducts();
        
//        double total_plus_reactive_counts = accumulate(treated_adducts.begin() + 1,
//                                                       treated_adducts.end(), 0.0);
        //double total_plus_counts = total_plus_reactive_counts + treated_adducts[0];
        vector<double> plus_channel_freq(treated_adducts.size(), 0);
//        for (size_t i = 1; i < treated_adducts.size(); i++)
//        {
//            plus_channel_freq[i] = (treated_adducts[i] / total_plus_reactive_counts);
//        }
        
//        double total_minus_reactive_counts = accumulate(untreated_adducts.begin() + 1,
//                                                        untreated_adducts.end(), 0.0);
        //double total_minus_counts = total_minus_reactive_counts + untreated_adducts[0];
        vector<double> minus_channel_freq(untreated_adducts.size(), 0);;
//        for (size_t i = 1; i < untreated_adducts.size(); i++)
//        {
//            minus_channel_freq[i] = (untreated_adducts[i] / total_minus_reactive_counts);
//        }

        assert (minus_channel_freq.size() == plus_channel_freq.size());
        
        // compute the betas
        for (size_t k = 0; k < plus_channel_freq.size(); k++)
        {
            double plus_denom = 0.0;
            double minus_denom = 0.0;
            
            // count the number of fragments that stopped RT up until reached
            // k (and include k).
            for (size_t i = k; i >=0 && i < plus_channel_freq.size(); i--)
            {
                plus_denom += treated_adducts[i];
                minus_denom += untreated_adducts[i];
            }
            
            if (plus_denom == 0 || minus_denom == 0)
            {
                continue;
            }
            
            double beta = treated_adducts[k]/plus_denom - untreated_adducts[k]/minus_denom; 
            
            if (untreated_adducts[k]/minus_denom != 1.0)
            {
                beta /= (1 - (untreated_adducts[k]/minus_denom));
            }
            else
            {
                beta = 0.0;
            }
            
            beta = max(0.0, beta);
            assert (!isinf(beta) && !isnan(beta));
            _betas[k] = beta;
        }
        
        double plus_denom = 0.0;
        double minus_denom = 0.0;
        for (size_t i = plus_channel_freq.size() - 1; i >=0 && i < plus_channel_freq.size(); i--)
        {
            plus_denom += treated_adducts[i];
            minus_denom += untreated_adducts[i];
        }
        
        //calculate c
        //double c = log(untreated_adducts[0] / minus_denom) - log(treated_adducts[0]/plus_denom);
        double c = 0.0;
        //JBL - does not matter that this sum runs backward according to the code (vs. paper) indexing convention
        for (size_t i = 1; i < plus_channel_freq.size(); ++i) 
        {
            if (_betas[i] > 0)
                c -= log(1 - _betas[i]);
        }
        
        if (isnan(c) || isinf(c))
        {
            // FIXME: drop to CE model?
            return;
        }
        
        //calculate theta
        for (size_t k = 1; k < plus_channel_freq.size(); k++)
        {
            double plus_denom = 0.0;
            double minus_denom = 0.0;
            
            for (size_t i = k; i >= 0 && i < plus_channel_freq.size(); i--)
            {
                plus_denom += treated_adducts[i];
                minus_denom += untreated_adducts[i];
            }
            
            if (plus_denom == 0 || minus_denom == 0)
            {
                continue;
            }
            
            double theta = log(1 - untreated_adducts[k]/minus_denom);
            theta -= log(1 - treated_adducts[k]/plus_denom);
            theta /= c;
            
            theta = max(0.0, theta);
            assert (!isinf(theta) && !isnan(theta));
            _thetas[k] = theta;
        }
        
        _c = c;
    }

    
    void print_adduct_counts(FILE* adducts_out)
    {
        const vector<int>& treated_adducts = _treated.adducts();
        const vector<int>& untreated_adducts = _untreated.adducts();
        
        vector<double> scaled_betas = _betas;
        double total_beta = accumulate(_betas.begin() + 1, _betas.end(), 0.0);
        for (size_t i = 1; i < scaled_betas.size(); ++i)
        {
            scaled_betas[i] /= total_beta;
        }
        
        //double total_theta = accumulate(_thetas.begin() + 1, _thetas.end(), 0.0); //JBL - not used
        total_beta = accumulate(scaled_betas.begin() + 1, scaled_betas.end(), 0.0);
        //fprintf(stderr, "%s: total theta = %lg\n", _name.c_str(), total_theta);
        //fprintf(stderr, "%s: total beta = %lg\n", _name.c_str(), total_beta);
        
        for (int i = 0; i < _profile_len; ++i)
        {
            if (i == 0)
            {
                fprintf(adducts_out, 
                        "%s\t%d\t*\t%d\t%d\t-\t-\t%lg\n", 
                        _name.c_str(),
                        i,
                        treated_adducts[i], 
                        untreated_adducts[i],
                        _c); 
            }
            else if (i < (int)_seq.length())
            {
                fprintf(adducts_out, 
                        "%s\t%d\t%c\t%d\t%d\t%lg\t%lg\t%lg\n", 
                        _name.c_str(),
                        i,
                        _seq[i-1],
                        treated_adducts[i], 
                        untreated_adducts[i],
                        scaled_betas[i],
                        _thetas[i],
                        _c); 
            }
        }
    }
    
private:
    string _name;
    string _seq;
    int _profile_len;
    //vector<int> _phys_cov;
    
    Adducts _treated;
    Adducts _untreated;
    vector<double> _thetas;
    vector<double> _betas;
    //vector<double> _normalized_thetas;
    double _c;
};

//wraps TargetProfile
//contains list of TargetProfiles, one at each priming spot

class RandomerTargetProfile
{
public:
    RandomerTargetProfile() {}
    RandomerTargetProfile(const string& name, const string& seq) :
        _name(name), 
        _seq(seq)
    {
        //constructing list of all possible starting positions for RT
        for (size_t i = 1; i <= seq.length(); ++i)
        {
            string subseq = seq.substr(0,i);
            _starts.push_back(TargetProfile(name, subseq));
        }
    }
    
    const string& name() const { return _name; }
    const string& seq() const { return _seq; }
    
    bool register_fragment(const MateHit& fragment, bool treated)
    {
        int rt_start_site = fragment.right(); // fragment is 0-indexed, adduct array is 1-indexed, so we don't actually need to subtract

        //if the rt_start_site is within the possible range, then register the fragment to that start site
        if (rt_start_site >= 1 && rt_start_site - 1 < (int)_starts.size())
        {                
            //JBL - not sure why the rt_start_site - 1 here
            //_starts[] is 0-indexed
            return _starts[rt_start_site - 1].register_fragment(fragment, treated); 
        }
        return false;
    }
    
    //calc reactivities profile for each start
    void calc_poisson_reactivities()
    {
        for (size_t i = 0; i < _starts.size(); ++i)
        {
            _starts[i].calc_poisson_reactivities();
        }
    }
    
    
    void print_adduct_counts(FILE* adducts_out)
    {
        //JBL - only print out all RT start site information if
        // all_RT_starts is set to true
        int j_start = _starts.size()-1;
        if (all_RT_starts) 
        {
            j_start = 0;
        }
        
        for (int j = j_start; j < (int)_starts.size(); ++j)
        {
            const TargetProfile& curr_start = _starts[j];
            
            const vector<int>& treated_adducts = curr_start.treated().adducts();
            const vector<int>& untreated_adducts = curr_start.untreated().adducts();
            
            vector<double> scaled_betas = curr_start.betas();
            double total_beta = accumulate(scaled_betas.begin() + 1, scaled_betas.end(), 0.0);
            if (total_beta == 0.0)
            {
                continue;
            }
            for (size_t i = 1; i < scaled_betas.size(); ++i)
            {
                scaled_betas[i] /= total_beta;
            }
            
            //double total_theta = accumulate(curr_start.thetas().begin() + 1, curr_start.thetas().end(), 0.0);
            //total_beta = accumulate(scaled_betas.begin() + 1, scaled_betas.end(), 0.0);
            //fprintf(stderr, "%s: total theta = %lg\n", _name.c_str(), total_theta);
            //fprintf(stderr, "%s: total beta = %lg\n", _name.c_str(), total_beta);
            
            for (int i = 0; i < j+1; ++i)
            {
                if (i == 0)
                {
                    fprintf(adducts_out, 
                            "%s\t%d\t%d\t*\t%d\t%d\t-\t-\t%lg\n", 
                            _name.c_str(),
                            j,
                            i,
                            treated_adducts[i], 
                            untreated_adducts[i],
                            curr_start.c()); 
                }
                else if (i < j)
                {
                    fprintf(adducts_out, 
                            "%s\t%d\t%d\t%c\t%d\t%d\t%lg\t%lg\t%lg\n", 
                            _name.c_str(),
                            j,
                            i,
                            _seq[i-1],
                            treated_adducts[i], 
                            untreated_adducts[i],
                            scaled_betas[i],
                            curr_start.thetas()[i],
                            curr_start.c()); 
                }
            }
        }
    }

    void print_consensus_reactivities(FILE* adducts_out)
    {
        for (int i = 0; i < (int)_consensus_thetas.size(); ++i)
        {
            if (i == 0)
            {
                fprintf(adducts_out,
                        "%s\t%lu\t%d\t*\t%d\t%d\t-\t-\t%lg\n",
                        _name.c_str(),
                        _seq.length(),
                        i,
                        _total_treated_adducts[i],
                        _total_untreated_adducts[i],
                        _consensus_c);
            }
            else
            {
                fprintf(adducts_out,
                        "%s\t%lu\t%d\t%c\t%d\t%d\t%lg\t%lg\t%lg\n",
                        _name.c_str(),
                        _seq.length(),
                        i,
                        _seq[i-1],
                        _total_treated_adducts[i],
                        _total_untreated_adducts[i],
                        _consensus_betas[i],
                        _consensus_thetas[i],
                        _consensus_c);
            }
        }
    }
    
    //JBL - Trapnell's averaging process over different start sites ocurrs here
    void calculate_consensus_reactivities()
    {
        calc_poisson_reactivities();
        
        _consensus_thetas = vector<double>(_seq.length() + 1, 0);
        _consensus_betas = vector<double>(_seq.length() + 1, 0);

        _total_treated_adducts = vector<int>(_seq.length() + 1, 0);
        _total_untreated_adducts = vector<int>(_seq.length() + 1, 0);
        
        _consensus_c = 0;
        
        //get SUM_n = Y_i,n, X_i,n
        for (int j = 0; j < (int)_starts.size(); ++j)
        {
            const TargetProfile& curr_start = _starts[j];
            
            const vector<int>& treated_adducts = curr_start.treated().adducts();
            const vector<int>& untreated_adducts = curr_start.untreated().adducts();
            
            for (size_t i = 0; i < treated_adducts.size(); ++i)
            {
                _total_treated_adducts[i] += treated_adducts[i];
                _total_untreated_adducts[i] += untreated_adducts[i];
            }
        }

        //performing weighted average here
        for (int j = 0; j < (int)_starts.size(); ++j)
        {
            const TargetProfile& curr_start = _starts[j];
            
            const vector<int>& treated_adducts = curr_start.treated().adducts();
            const vector<int>& untreated_adducts = curr_start.untreated().adducts();
            
            const vector<double>& curr_thetas = curr_start.thetas();
            const vector<double>& curr_betas = curr_start.betas();

            for (size_t i = 0; i < curr_thetas.size(); ++i)
            {
                double total_site_frags = (_total_treated_adducts[i] + _total_untreated_adducts[i]);
                double curr_site_frags = (treated_adducts[i] + untreated_adducts[i]);
                if (total_site_frags > 0)
                {
                    _consensus_thetas[i] += curr_thetas[i] * (curr_site_frags / total_site_frags);
                    _consensus_betas[i] += curr_betas[i] * (curr_site_frags / total_site_frags);
                }
                else
                {
                    _consensus_thetas[i] = 0;
                    _consensus_betas[i] = 0;
                }
            }
        }
        
        double total_beta = accumulate(_consensus_betas.begin() + 1, _consensus_betas.end(), 0.0);
        if (total_beta > 0.0)
        {
            for (size_t i = 1; i < _consensus_betas.size(); ++i)
            {
                _consensus_betas[i] /= total_beta;
            }
        }
        
        double total_treated_adducts_across_sites = accumulate(_total_treated_adducts.begin() + 1, _total_treated_adducts.end(), 0.0);
        for (int j = 0; j < (int)_starts.size(); ++j)
        {
            const TargetProfile& curr_start = _starts[j];
            if (curr_start.c() > 0 && total_treated_adducts_across_sites > 0)
            {
                const vector<int>& treated_adducts = curr_start.treated().adducts();
                double total_treated_adducts_in_curr_site = accumulate(treated_adducts.begin() + 1, treated_adducts.end(), 0.0);
                _consensus_c += curr_start.c() * total_treated_adducts_in_curr_site / total_treated_adducts_across_sites;
            }
        }
    }
    
    
private:
    string _name;
    string _seq;
    
    vector<TargetProfile> _starts;
    
    vector<double> _consensus_thetas;
    vector<double> _consensus_betas;
    
    vector<int> _total_treated_adducts;
    vector<int> _total_untreated_adducts;

    double _consensus_c;
    
};

map<RefID, RandomerTargetProfile> targets_by_id;

void register_targets(RefSequenceTable rt,
                      FILE* target_fasta)
{
    while(!feof(target_fasta))
    {
        Read curr_target;
        if (!next_fasta_record(target_fasta, 
                               curr_target.name, 
                               curr_target.seq))
        {
            break;
        }
        
        string name = curr_target.name;
        RefID refid = rt.get_id(name, NULL);
      
        map<RefID,RandomerTargetProfile>::iterator itr = targets_by_id.find(refid);
        if (itr == targets_by_id.end())
        {
            targets_by_id[refid] = RandomerTargetProfile(name, curr_target.seq);
            assert(!(targets_by_id[refid].name().empty()));
        }
        else
        {
            assert (itr->second.seq().length() == curr_target.seq.length()); 
        }

    }
}

void process_fragments(FragmentFactory& fragment_factory, 
                       map<RefID, RandomerTargetProfile>& targets_by_id,
                       int& processed_fragments,
                       int& num_frags,
                       int& num_kept_frags,
                       bool treated)
{
    MateHit fragment;
    while(fragment_factory.next_fragment(fragment))
    {
        processed_fragments++;
        map<RefID, RandomerTargetProfile>::iterator itr;
        if (fragment.ref_id() == 0)
            return;
        itr = targets_by_id.find(fragment.ref_id());
        if (itr == targets_by_id.end())
        {
            fprintf(stderr, "warning: unknown target Ref ID encountered in Bowtie map\n");
            //continue;
            exit(1);
        }
        
        RandomerTargetProfile& target = itr->second;
        bool kept = target.register_fragment(fragment, treated);
        
        num_frags++;
        if (kept)
            num_kept_frags++;
		
        //int length = fragment.right() - fragment.left(); //JBL - unused
//        if (length > 0)
//        {
//			if (treated)
//			{				
//				if (treated_frag_length_dist.size() < length + 1)
//					treated_frag_length_dist.resize(length + 1);
//				treated_frag_length_dist[length]++;
//			}
//			else
//			{				
//				if (untreated_frag_length_dist.size() < length + 1)
//					untreated_frag_length_dist.resize(length + 1);
//				untreated_frag_length_dist[length]++;
//			}
//        }
    }
    
}

void driver(FILE* target_fasta, FILE* treated_sam_hits_file, FILE* untreated_sam_hits_file)
{
    RefSequenceTable rt(true, false);

    register_targets(rt, target_fasta);
    
    ReadTable it;
    SAMHitFactory hs(it, rt);
    FragmentFactory treated_frag_factory(hs, treated_sam_hits_file);
    FragmentFactory untreated_frag_factory(hs, untreated_sam_hits_file);
    
    int processed_fragments = 0;

    vector<int> treated_frag_length_dist;
	vector<int> untreated_frag_length_dist;
    
    int num_kept_treated_frags = 0;
    int num_kept_untreated_frags = 0;
	
	int num_treated_frags = 0;
    int num_untreated_frags = 0;
	
    process_fragments(treated_frag_factory, 
                      targets_by_id,
                      processed_fragments,
                      num_treated_frags,
                      num_kept_treated_frags,
                      true);
    
    process_fragments(untreated_frag_factory, 
                      targets_by_id,
                      processed_fragments,
                      num_untreated_frags,
                      num_kept_untreated_frags,
                      false);
    
    string adducts_out_name = output_dir + "/" + "reactivities.out";
    FILE* adducts_out = fopen(adducts_out_name.c_str(), "w");
    if (adducts_out == NULL)
    {
        fprintf(stderr, "Error: could not open %s for writing\n", 
                adducts_out_name.c_str());
        exit(1);
    }
    
    //printing reactivities for each target
    fprintf(adducts_out, "sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n");
    map<RefID, RandomerTargetProfile>::iterator itr;
    for (itr = targets_by_id.begin(); itr != targets_by_id.end(); ++itr)
    {
        RandomerTargetProfile target = itr->second;
        //string target_prefix = output_dir + "/" + target.name();
        //string adducts_out_name = target_prefix + ".adducts";
        
        if (compute_consensus_reactivities)
        {
            target.calculate_consensus_reactivities();
            target.print_consensus_reactivities(adducts_out);
        }else{
            target.calc_poisson_reactivities();
            target.print_adduct_counts(adducts_out);
        }
    }
    fclose(adducts_out);
    
    //JBL START - not sure this code block is correctly implemented/relevant now
    string treated_library_file = output_dir + "/treated_library_length.hist";
    FILE* treated_library_out = fopen(treated_library_file.c_str(), "w");
    if (treated_library_out == NULL)
    {
        fprintf(stderr, "Error: could not open %s for writing\n", 
                treated_library_file.c_str());
        exit(1);
    }
	
	string untreated_library_file = output_dir + "/untreated_library_length.hist";
    FILE* untreated_library_out = fopen(untreated_library_file.c_str(), "w");
    if (untreated_library_out == NULL)
    {
        fprintf(stderr, "Error: could not open %s for writing\n", 
                untreated_library_file.c_str());
        exit(1);
    }
    
    
    string stats_file = output_dir + "/mapping_stats.txt";
    FILE* stats_out = fopen(stats_file.c_str(), "w");
    if (stats_out == NULL)
    {
        fprintf(stderr, "Error: could not open %s for writing\n", 
                stats_file.c_str());
        exit(1);
    }
    
    fprintf(treated_library_out, "length\tnum_frags\n");
    for (int i = 0; i < (int)treated_frag_length_dist.size(); ++i)
    {
        fprintf(treated_library_out, "%d\t%d\n", i, treated_frag_length_dist[i]); 
    }
	
	fprintf(untreated_library_out, "length\tnum_frags\n");
    for (int i = 0; i < (int)untreated_frag_length_dist.size(); ++i)
    {
        fprintf(untreated_library_out, "%d\t%d\n", i, untreated_frag_length_dist[i]); 
    }
    //JBL END
    
//    fprintf(stats_out, "Processed %d properly paired fragments, kept %d\n", 
//            processed_fragments, num_kept_frags);
    fprintf(stderr, "Processed %d properly paired fragments, kept %d/%d (%f %%) treated, %d/%d (%f %%)  untreated\n", 
			processed_fragments, 
			num_kept_treated_frags,
			num_treated_frags,
			(float)num_kept_treated_frags/(float)num_treated_frags,
			num_kept_untreated_frags, 
			num_untreated_frags,
			(float)num_kept_untreated_frags/(float)num_untreated_frags);
    
//    fprintf(stats_out, "target\ttreated_fragments\tuntreated_fragments\ttreated_stops\tuntreated_stops\n");
//    for (map<RefID, TargetProfile>::iterator itr = targets_by_id.begin(); 
//         itr != targets_by_id.end(); ++itr) 
//    {
//        fprintf(stats_out, "%s\t%d\t%d\t%d\t%d\n", 
//                itr->second.name().c_str(),
//                itr->second.treated().total_fragments(),
//                itr->second.untreated().total_fragments(),
//                itr->second.treated().total_adducts(),
//                itr->second.untreated().total_adducts());
//    }
    
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv, print_usage);
    if (parse_ret)
        return parse_ret;
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string targets_file_name = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string treated_sam_hits_file_name = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string untreated_sam_hits_file_name = argv[optind++];
	
    // Open the approppriate files
    FILE* fasta_targets_file = fopen(targets_file_name.c_str(), "r");
    if (fasta_targets_file == NULL)
    {
        fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
                targets_file_name.c_str());
        exit(1);
    }
    
    FILE* treated_sam_hits_file = fopen(treated_sam_hits_file_name.c_str(), "r");
    if (treated_sam_hits_file == NULL)
    {
        fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
                treated_sam_hits_file_name.c_str());
        exit(1);
    }
    
    FILE* untreated_sam_hits_file = fopen(untreated_sam_hits_file_name.c_str(), "r");
    if (untreated_sam_hits_file == NULL)
    {
        fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
                untreated_sam_hits_file_name.c_str());
        exit(1);
    }
	
	srand48(time(NULL));
	
    driver(fasta_targets_file, treated_sam_hits_file, untreated_sam_hits_file);
	
	return 0;
}

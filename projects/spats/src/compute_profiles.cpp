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
        if (adduct_site >= 0 && adduct_site < adduct_counts.size())
        {
            int end = min((int)adduct_counts.size(), fragment.right());
            if (end == adduct_counts.size() - 1)
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

struct TargetProfile
{
    TargetProfile() {}
    
    TargetProfile(const string& name, const string& seq) : 
        _name(name), 
        _seq(seq),
        _profile_len(_seq.length() + 1),
        _treated(_profile_len),
        _untreated(_profile_len),
        _thetas(_seq.length()),
        _normalized_thetas(_seq.length()){}
    
    const string& name() const { return _name; }
    const string& seq() const { return _seq; }
    
    const Adducts treated() const { return _treated; }
    const Adducts untreated() const { return _untreated; }
    
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
        if (_profile_len == 0)
            return;
        
        if (_name == "target_WT")
        {
            int a =4;
        }
        
        const vector<int>& treated_adducts = _treated.adducts();
        const vector<int>& untreated_adducts = _untreated.adducts();
        
        double X_0 = treated_adducts[0];
        double total_plus_reactive_counts = accumulate(treated_adducts.begin() + 1,
                                                       treated_adducts.end(), 0.0);
        double total_plus_counts = total_plus_reactive_counts + treated_adducts[0];
        vector<double> plus_channel_freq;
        for (size_t i = 1; i < treated_adducts.size(); i++)
        {
            plus_channel_freq.push_back(treated_adducts[i] / total_plus_reactive_counts);
        }
        
        double total_minus_reactive_counts = accumulate(untreated_adducts.begin() + 1,
                                                        untreated_adducts.end(), 0.0);
        double total_minus_counts = total_minus_reactive_counts + untreated_adducts[0];
        vector<double> minus_channel_freq;
        for (size_t i = 1; i < untreated_adducts.size(); i++)
        {
            minus_channel_freq.push_back(untreated_adducts[i] / total_minus_counts);
        }
        
        if (total_minus_counts <= 0.0)
        {
            fprintf(stderr, "Error: no full length fragments in the untreated channel for %s\n", _name.c_str());
            return;
        }
        
        double p_0_hat = untreated_adducts[0] / (total_minus_counts);
        
        double scaled_X_0 = X_0 / p_0_hat;
        
        double cap_C_estimate = (total_plus_counts) / scaled_X_0;
        
        //assert (cap_C_estimate > 0.0);
        
        double c_estimate = log(cap_C_estimate);
        
        assert (minus_channel_freq.size() == plus_channel_freq.size());
        
        double K = p_0_hat / (cap_C_estimate - p_0_hat);
        
        for (int i = 0; i < (int)plus_channel_freq.size(); ++i)
        {
            double P = 0.0;
            double M = 0.0;
            for (int j = i - 1; j >= 0; --j)
            {
                P += (plus_channel_freq[j]);
                M += (minus_channel_freq[j]);
            }
            double p = 0.0;
            if (P + K > 0)
                p = (plus_channel_freq[i] / (P + K));
            
            double m = 0.0;
            if (M + p_0_hat > 0)
                m = (minus_channel_freq[i] / (M + p_0_hat));
            
            double log_p = log(1.0 + p);
            double log_m = log(1.0 + m);
            double theta = (1.0 / c_estimate) * (log_p - log_m);
            _thetas[i] = theta;
        }
        
        double delta = 1.0;
        
        for (int i = 0; i < _thetas.size(); ++i)
        {
            if (_thetas[i] < 0)
                delta -= _thetas[i];
        }
        
        for (int i = 0; i < _thetas.size(); ++i)
        {
            if (_thetas[i] < 0)
            {
                _normalized_thetas[i] = 0;
            }
            else 
            {
                _normalized_thetas[i] = _thetas[i] / delta;
            }
        }
    }
    
    void print_adduct_counts(FILE* adducts_out)
    {
        fprintf(adducts_out, "five_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\n");
        
        const vector<int>& treated_adducts = _treated.adducts();
        const vector<int>& untreated_adducts = _untreated.adducts();
        
        for (int i = 0; i < _profile_len; ++i)
        {
            if (i == 0)
            {
                fprintf(adducts_out, 
                        "%d\t*\t%d\t%d\t-\t-\n", 
                        i,
                        treated_adducts[i], 
                        untreated_adducts[i]); 
            }
            else if (i < _seq.length())
            {
                fprintf(adducts_out, 
                        "%d\t%c\t%d\t%d\t%lg\t%lg\n", 
                        i,
                        _seq[i-1],
                        treated_adducts[i], 
                        untreated_adducts[i],
                        _thetas[i-1],
                        _normalized_thetas[i-1]); 
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
    vector<double> _normalized_thetas;
};

map<RefID, TargetProfile> targets_by_id;

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
      
        map<RefID,TargetProfile>::iterator itr = targets_by_id.find(refid);
        if (itr == targets_by_id.end())
        {
            targets_by_id[refid] = TargetProfile(name, curr_target.seq);
            assert(!(targets_by_id[refid].name().empty()));
        }
        else
        {
            assert (itr->second.seq().length() == curr_target.seq.length()); 
        }

    }
}

void process_fragments(FragmentFactory& fragment_factory, 
                       map<RefID, TargetProfile>& targets_by_id,
                       int& processed_fragments,
                       int& num_frags,
                       int& num_kept_frags,
                       bool treated)
{
    MateHit fragment;
    while(fragment_factory.next_fragment(fragment))
    {
        processed_fragments++;
        map<RefID, TargetProfile>::iterator itr;
        if (fragment.ref_id() == 0)
            return;
        itr = targets_by_id.find(fragment.ref_id());
        if (itr == targets_by_id.end())
        {
            fprintf(stderr, "Error: unknown target Ref ID encountered in Bowtie map\n");
            exit(1);
        }
        
        TargetProfile& target = itr->second;
        bool kept = target.register_fragment(fragment, treated);
        
        num_frags++;
        if (kept)
            num_kept_frags++;
		
        int length = fragment.right() - fragment.left();
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
    
    map<RefID, TargetProfile>::iterator itr;
    for (itr = targets_by_id.begin(); itr != targets_by_id.end(); ++itr)
    {
        TargetProfile target = itr->second;
        string target_prefix = output_dir + "/" + target.name();
        string adducts_out_name = target_prefix + ".adducts";
        FILE* adducts_out = fopen(adducts_out_name.c_str(), "w");
        if (adducts_out == NULL)
        {
            fprintf(stderr, "Error: could not open %s for writing\n", 
                    adducts_out_name.c_str());
            exit(1);
        }
        target.calc_poisson_reactivities();
        target.print_adduct_counts(adducts_out);
        fclose(adducts_out);
    }
    
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
    
    
//    fprintf(stats_out, "Processed %d properly paired fragments, kept %d\n", 
//            processed_fragments, num_kept_frags);
    fprintf(stderr, "Processed %d properly paired fragments, kept %d/%d ( %f\\% ) treated, %d/%d ( %f\\% ) untreated\n", 
			processed_fragments, 
			num_kept_treated_frags,
			num_treated_frags,
			(float)num_kept_treated_frags/(float)num_treated_frags,
			num_kept_untreated_frags, 
			num_untreated_frags,
			(float)num_kept_untreated_frags/(float)num_untreated_frags);
    
    fprintf(stats_out, "target\ttreated_fragments\tuntreated_fragments\ttreated_stops\tuntreated_stops\n");
    for (map<RefID, TargetProfile>::iterator itr = targets_by_id.begin(); 
         itr != targets_by_id.end(); ++itr) 
    {
        fprintf(stats_out, "%s\t%d\t%d\t%d\t%d\n", 
                itr->second.name().c_str(),
                itr->second.treated().total_fragments(),
                itr->second.untreated().total_fragments(),
                itr->second.treated().total_adducts(),
                itr->second.untreated().total_adducts());
    }
    
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
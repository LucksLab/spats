import os

from db import PairDB
from diagram import diagram
from pair import Pair
from parse import reactivities_parse
from spats import Spats


def compare_v102(compare_results_path, targets_path, treated_mask, untreated_mask, data_r1_path, data_r2_path, spats_out_path, beta_theta_compare_threshold = 0.001):

    db_path = os.path.join(compare_results_path, "pairs.db")
    db = PairDB(db_path)
    db.wipe()
    count = db.parse(data_r1_path, data_r2_path)
    db.index()
    print("Parsed and indexed {} data records.".format(count))

    db.add_v102_comparison(targets_path, spats_out_path)

    spats = Spats()
    spats.run.minimum_target_match_length = 10
    spats.run.masks = [ treated_mask, untreated_mask ]
    spats.addTargets(targets_path)
    print("Current build processing pair data...")
    spats.writeback_results = True
    spats.process_pair_db(db)

    all_profiles = spats.compute_profiles()
    spats.write_reactivities(os.path.join(compare_results_path, 'rx.out'))

    print("current missing from v102: {}".format(len(db.our_pairs_missing_from_v102())))
    print("v102 missing from current: {}".format(len(db.v102_pairs_missing_from_ours())))
    print("Ours differing from v102: {}".format(len(db.our_pairs_differing_from_v102())))

    print("Reasons for v102 missing from current:")
    spats.run.debug = True
    pair = Pair()
    diagrams = 0
    reasons = 0
    for reason, count in db.v102_pairs_missing_from_ours_reasons():
        reasons += 1
        print("  {}: {}".format(count, reason))
        with open(os.path.join(compare_results_path, 'v102_delta_{}.out'.format(reason)), 'wb') as outfile:
            spats.run.log = outfile
            for p in db.our_pairs_nonmatching_v102_for_reason(reason):
                outfile.write("="*100 + "\n")
                outfile.write("pair rowid: {}, v102 numeric id: {}, site: {}, NOMASK_R1: {}, NOMASK_R2: {}\n".format(p[0], p[4], p[5], p[6], p[7]))
                pair.set_from_data(str(p[1]), str(p[2]), str(p[3]))
                spats.process_pair(pair)
                outfile.write("-"*100 + "\n")
                outfile.write(diagram(pair) + "\n")
                diagrams += 1
    print("Wrote {} differing diagrams (for {} reasons) to {}".format(diagrams, reasons, compare_results_path))

    sites_diff = 0
    with open(os.path.join(compare_results_path, 'reactivities_diff.out'), 'wb') as outfile:
        outfile.write("Differences > {} in beta/theta (current first / v102 second)\n".format(beta_theta_compare_threshold))
        outfile.write("         beta     theta     treated   untreated\n")
        only_01_diffs = ""
        for entry in reactivities_parse(os.path.join(spats_out_path, 'reactivities.out')):
            # return list of (target, rt_start, site, nuc, treated_count, untreated_count, beta, theta, c)
            site = int(entry[2])
            if 0 == site:
                continue
            target_name = str(entry[0])
            profiles = all_profiles.profilesForTargetNamed(target_name)
            beta = float(entry[6])
            theta = float(entry[7])
            tc = int(entry[4])
            uc = int(entry[5])
            our_beta = profiles.betas[site]
            our_theta = profiles.thetas[site]
            our_tc = profiles.treated_counts[site]
            our_uc = profiles.untreated_counts[site]
            if abs(beta - our_beta) > beta_theta_compare_threshold or abs(theta - our_theta) > beta_theta_compare_threshold:
                report =  "\ncur    {:.6f}  {:.6f}    {}         {}    \t  {} :: {}".format(our_beta, our_theta, our_tc, our_uc, target_name, site)
                report += "\nv102   {:.6f}  {:.6f}    {}         {}\n".format(beta, theta, tc, uc)
                if set([tc, uc, our_tc, our_uc]) <= set([0, 1]):
                    only_01_diffs += report
                else:
                    outfile.write(report)
                sites_diff += 1
        if only_01_diffs:
            outfile.write("\n" + ("-" * 50) + "\n")
            outfile.write(only_01_diffs)
    if sites_diff > 0:
        print("Wrote {} sites with beta/theta diff > {} to {}".format(sites_diff, beta_theta_compare_threshold, compare_results_path))
    else:
        print("All reactivities within {} threshold".format(beta_theta_compare_threshold))

def diagram_case(db_path, rowid, targets_path, treated_mask, untreated_mask):
    db = PairDB(db_path)
    pair = Pair()
    p = db.diagram_info(rowid)
    print("pair rowid: {}, v102 numeric id: {}, site: {}, NOMASK_R1: {}, NOMASK_R2: {}".format(p[0], p[4], p[5], p[6], p[7]))
    pair.set_from_data(str(p[1]), str(p[2]), str(p[3]))

    spats = Spats()
    spats.run.minimum_target_match_length = 10
    spats.run.masks = [ treated_mask, untreated_mask ]
    spats.addTargets(targets_path)
    spats.process_pair(pair)
    print(diagram(pair))
    

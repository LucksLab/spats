
import os
import random
import sqlite3
import sys

from parse import FastFastqParser

class PairDB(object):

    def __init__(self, database_path = None):
        self._dbpath = database_path or ":memory:"
        self.conn = sqlite3.connect(self._dbpath)
        self.show_progress_every = False

    def wipe(self):
        self.conn.execute("DROP TABLE IF EXISTS pair")
        self.conn.execute("DROP INDEX IF EXISTS r1_idx")
        self.conn.execute("DROP INDEX IF EXISTS r2_idx")
        self.conn.execute("DROP INDEX IF EXISTS pair_id_idx")
        self.conn.execute("DROP TABLE IF EXISTS tag")
        self.conn.execute("DROP TABLE IF EXISTS result")
        self.conn.execute("DROP TABLE IF EXISTS result_set")
        self._wipe_v102()

    def _create(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS pair (r1 TEXT, r2 TEXT, identifier TEXT)")

    def load_and_index(self, target_path, r1_path, r2_path):
        self.add_targets_table(target_path)
        self.parse(r1_path, r2_path)
        self._cache_unique_pairs()


    # we might be able to speed this up by:
    #  - some sqlite tricks like: conn.execute("PRAGMA synchronous=OFF"), conn.execute("PRAGMA cache_size=16000"), etc.
    #  - parsing in one thread and inserting in another
    # but for now this seems fast enough (~300k pairs/s?)
    def parse(self, r1_path, r2_path):
        conn = self.conn
        conn.execute("DROP INDEX IF EXISTS r1_idx")
        conn.execute("DROP INDEX IF EXISTS r2_idx")
        self._create()
        total = 0
        show_progress_every = self.show_progress_every
        next_report_count = show_progress_every
        with FastFastqParser(r1_path, r2_path) as parser:
            while True:
                pairs, count = parser.read(16384)
                total += count
                if not pairs:
                    break
                conn.executemany("INSERT INTO pair (identifier, r1, r2) VALUES (?, ?, ?)", pairs)
                if show_progress_every:
                    next_report_count -= count
                    if next_report_count < 0:
                        next_report_count = show_progress_every
                        sys.stdout.write('.')
                        sys.stdout.flush()
        conn.commit()
        return total

    def index(self):
        self.conn.commit()
        self.conn.execute("CREATE INDEX IF NOT EXISTS r1_idx ON pair (r1)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS r2_idx ON pair (r2)")

    def _fetch_one(self, query, args = []):
        res = self.conn.execute(query, args).fetchone()
        return res[0] if res else None

    def count(self):
        return self._fetch_one("SELECT COUNT(*) from pair")

    def unique_r1(self):
        self.index()
        return self._fetch_one("SELECT COUNT(distinct r1) from pair")

    def unique_r2(self):
        self.index()
        return self._fetch_one("SELECT COUNT(distinct r2) from pair")

    def _cache_unique_pairs(self):
        self.index()
        self.conn.execute("CREATE TABLE IF NOT EXISTS unique_pair (pair_id INT, multiplicity INT)")
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS unique_pair_idx ON unique_pair (pair_id)")
        num_already = self._fetch_one("SELECT 1 FROM unique_pair limit 1")
        if not num_already:
            self.conn.execute("INSERT INTO unique_pair (pair_id, multiplicity) SELECT rowid, COUNT(rowid) FROM pair GROUP BY r1||r2")
            self.conn.commit()

    def pair_length(self):
        return len(self._fetch_one("SELECT r1 FROM pair LIMIT 1"))

    def unique_pairs(self):
        self._cache_unique_pairs()
        return self._fetch_one("SELECT COUNT(*) from unique_pair")

    def max_r1(self):
        self.index()
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r1 order by cnt desc limit 1")

    def max_r2(self):
        self.index()
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r2 order by cnt desc limit 1")

    def _batch_results(self, query, batch_size, args = []):
        offset = 0
        while True:
            use_query = query.format(offset, batch_size)
            batch = sqlite3.connect(self._dbpath).execute(use_query, args).fetchall()
            if batch:
                offset = (int(batch[-1][4]) + 1)
                yield batch
            else:
                return

    def unique_pairs_with_counts(self, batch_size = 0):
        self._cache_unique_pairs()
        return self._batch_results('''SELECT u.multiplicity, p.r1, p.r2, p.rowid, u.rowid
                                      FROM unique_pair u JOIN pair p ON p.rowid=u.pair_id
                                      WHERE u.rowid >= {} ORDER BY u.rowid LIMIT {}''', batch_size)

    def unique_pairs_with_counts_and_no_results(self, result_set_id, batch_size = 0):
        self._cache_unique_pairs()
        self.index_results()
        return self._batch_results('''SELECT u.multiplicity, p.r1, p.r2, p.rowid, u.rowid
                                      FROM unique_pair u
                                      JOIN pair p ON p.rowid=u.pair_id
                                      LEFT JOIN result r ON r.set_id=? AND r.pair_id=u.pair_id
                                      WHERE r.rowid IS NULL AND u.rowid >= {}
                                      ORDER BY u.rowid LIMIT {}''', batch_size, (result_set_id,))

    def all_pairs(self, batch_size = 0):
        return self._batch_results("SELECT 1, r1, r2, rowid FROM pair WHERE rowid >= {} ORDER BY rowid LIMIT {}", batch_size)

    def add_targets_table(self, targets_path):
        from parse import fasta_parse
        self.conn.execute("DROP TABLE IF EXISTS target")
        self.conn.execute("DROP INDEX IF EXISTS target_idx")
        self.conn.execute("CREATE TABLE target (name TEXT, seq TEXT)")
        self.conn.execute("CREATE INDEX target_idx ON target (name)")
        target_list = fasta_parse(targets_path)
        self.conn.executemany('INSERT INTO target (name, seq) VALUES (? , ?)', target_list)
        return { name : seq for name, seq in target_list }

    def targets(self):
        return [ (r[0], r[1], r[2]) for r in self.conn.execute("SELECT name, seq, rowid FROM target") ]

    # results storing
    def _prepare_results(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS result (set_id INT, pair_id INT, target INT, site INT, mask TEXT, multiplicity INT, failure TEXT)")
        self.conn.execute("DROP INDEX IF EXISTS pair_result_idx")

    def add_result_set(self, set_name, resume_processing = False):
        self._prepare_results()
        self.conn.execute("CREATE TABLE IF NOT EXISTS result_set (set_id TEXT)")
        rid = self.result_set_id_for_name(set_name)
        if rid is not None and not resume_processing:
            self.conn.execute("DELETE FROM result WHERE set_id=?", (rid,))
            self.conn.execute("DELETE FROM result_set WHERE rowid=?", (rid,))
            rid = None
        if rid is None:
            self.conn.execute("INSERT INTO result_set (set_id) VALUES (?)", (set_name,))
        self.conn.commit()
        return self.result_set_id_for_name(set_name)

    def result_set_id_for_name(self, set_name):
        return self._fetch_one("SELECT rowid FROM result_set WHERE set_id=?", (set_name,))

    # results a list of (rowid, target_name, site, mask, multiplicity, failure, [optional: tags])
    def add_results(self, result_set_id, results):
        # grab a new connection since this might be in a new process (due to multiprocessing)
        #print " --> Thd in AR {}".format(self.worker_id)
        has_tags = (len(results[0]) > 6)
        conn = sqlite3.connect(self._dbpath)
        stmt = '''INSERT INTO result (set_id, pair_id, target, site, mask, multiplicity, failure)
                  VALUES ({}, ?, ?, ?, ?, ?, ?)'''.format(result_set_id)
        cursor = conn.executemany(stmt, [ r[0:6] for r in results ] if has_tags else results)
        if cursor.rowcount != len(results):
            raise Exception("some results failed to add: {} / {}".format(cursor.rowcount, len(results)))
        if has_tags:
            rstmt_template = 'INSERT INTO result_tag (result_id, tag_id) VALUES ({}, ?)'
            for res in results:
                rid = conn.execute("SELECT rowid FROM result WHERE set_id=? AND pair_id=?", (result_set_id, res[0])).fetchone()[0]
                conn.executemany(rstmt_template.format(rid), [ (t,) for t in res[6] ])
        conn.commit()

    def index_results(self):
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS pair_result_idx ON result (set_id, pair_id)")

    def num_results(self, result_set_name):
        return self._fetch_one("SELECT count(*) FROM result r JOIN result_set n ON r.set_id = n.rowid WHERE n.set_id = ?", (result_set_name,))

    def differing_results(self, result_set_name_1, result_set_name_2):
        self.index_results()
        id1 = self.result_set_id_for_name(result_set_name_1)
        id2 = self.result_set_id_for_name(result_set_name_2)
        return self.conn.execute('''SELECT p.rowid, p.r1, p.r2, s1.target, s1.site, s1.failure, s2.target, s2.site, s2.failure
                                    FROM result s1 JOIN result s2 ON s2.set_id=? AND s2.pair_id=s1.pair_id
                                    JOIN pair p ON p.rowid=s1.pair_id
                                    WHERE s1.set_id=? AND (s1.site != s2.site OR s1.target != s2.target)''', (id2, id1))

    # results analysis
    def setup_tags(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS tag (name TEXT)")
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS tag_name_idx ON tag (name)")
        self.conn.execute("CREATE TABLE IF NOT EXISTS result_tag (result_id INT, tag_id INT)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS result_tag_idx ON result_tag (tag_id)")

    def add_tags(self, tags):
        existing = self.tagmap()
        for t in tags:
            if -1 == existing.get(t, -1):
                self.conn.execute("INSERT INTO tag (name) VALUES (?)", (t,) )
        self.conn.commit()

    def tagmap(self):
        return { str(res[1]) : int(res[0]) for res in self.conn.execute("SELECT rowid, name FROM tag") }

    def count_tags(self):
        self.conn.execute("DROP TABLE IF EXISTS tag_count")
        self.conn.execute("CREATE TABLE tag_count (tag_id INT, count INT)")
        self.conn.execute('''INSERT INTO tag_count
                             SELECT t.rowid, SUM(r.multiplicity)
                             FROM result_tag rt
                             JOIN tag t ON t.rowid=rt.tag_id
                             JOIN result r ON r.rowid=rt.result_id''')


    #v102 delta analysis
    def _create_v102(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS v102 (pair_id INT, numeric_id INT, target INT, site INT, nomask_r1 TEXT, nomask_r2 TEXT, mask TEXT)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS pair_id_idx ON pair (identifier)")

    def _wipe_v102(self):
        self.conn.execute("DROP TABLE IF EXISTS v102")
        self.conn.execute("DROP INDEX IF EXISTS v102_numeric_idx")
        self.conn.execute("DROP INDEX IF EXISTS pair_id_idx")
        self.conn.execute("DROP INDEX IF EXISTS v102_site_idx")

    def add_v102_comparison(self, targets_path, output_path):
        import sys
        from parse import reactivities_parse, SamParser
        self.index()
        self._wipe_v102()
        self._create_v102()
        total = 0
        conn = self.conn

        # create the targets table
        targets = self.add_targets_table(targets_path)

        # parse NOMASK to get the ID<->numeric_id map
        with FastFastqParser(os.path.join(output_path, 'NOMASK_1.fq'), os.path.join(output_path, 'NOMASK_2.fq')) as parser:
            while True:
                pairs, count = parser.read_nomask(16384)
                total += count
                if not pairs:
                    break
                conn.executemany('''INSERT INTO v102
                                           (pair_id, numeric_id, target, site, nomask_r1, nomask_r2, mask)
                                    SELECT  pair.rowid,       ?,   NULL,   -1,         ?,         ?, NULL
                                    FROM pair WHERE pair.identifier=?''', pairs)
                sys.stdout.write('.')
                sys.stdout.flush()
        conn.commit()
        conn.execute("CREATE INDEX IF NOT EXISTS v102_numeric_idx ON v102 (numeric_id)")
        check = self._fetch_one("SELECT COUNT(*) from v102")
        if check != total:
            raise Exception("some v102 pairs did not match? {} != {}".format(check, total))
        print "\nParsed {} records from NOMASK_{{1,2}}.fq".format(total)

        # now, parse sam to get the numeric_id->site map
        total = 0
        mask_totals = {}
        masks = ('RRRY', 'YYYR')
        for mask in masks:
            mask_total = 0
            with SamParser(os.path.join(output_path, mask + '.sam'), targets) as parser:
                while True:
                    pairs, count = parser.read(16384, mask)
                    if not pairs:
                        break
                    total += count
                    mask_total += count
                    conn.executemany('UPDATE v102 SET target=(SELECT rowid FROM target WHERE name=?), site=?, mask=? WHERE numeric_id=?', pairs)
                    sys.stdout.write('.')
                    sys.stdout.flush()
            mask_totals[mask] = mask_total
        conn.commit()
        conn.execute("CREATE INDEX IF NOT EXISTS v102_site_idx ON v102 (site)")
        print "\nAdded {} records from *.sam:".format(total)
        for mask in masks:
            print "  {}.sam: {}".format(mask, mask_totals[mask])

        # now, compare against reactivities.out to verify our counts
        all_counts = { "{}::{}::{}".format(r[0],r[1],r[2]) : r[3] for r in self.conn.execute('''SELECT t.name, v.mask, v.site, count(v.rowid)
                                                                                                FROM v102 v
                                                                                                JOIN target t ON v.target=t.rowid
                                                                                                GROUP BY t.name, v.site, v.mask''') }
        count = 0
        for entry in reactivities_parse(os.path.join(output_path, 'reactivities.out')):
            # (target, rt_start, site, nuc, treated_count, untreated_count, beta, theta, c)
            key = "{}::{}::{}".format(entry[0], 'RRRY', int(entry[2]))
            if all_counts.get(key, 0) != int(entry[4]):
                raise Exception("Treated db does not match reactivities: key={}, db={}, r={}".format(key, all_counts.get(key, 0), int(entry[4])))
            key = "{}::{}::{}".format(entry[0], 'YYYR', int(entry[2]))
            if all_counts.get(key, 0) != int(entry[5]):
                raise Exception("Untreated db does not match reactivities: key={}, db={}, r={}".format(key, all_counts.get(key, 0), int(entry[5])))
            count += 1
        if not count:
            raise Exception("no sites found in reactivities.out?")
        print "reactivities.out check pass ({} sites).".format(count)
        self.add_v102_to_results()

    def add_v102_to_results(self):
        v102_id = self.add_result_set("v102")
        self.conn.execute('''INSERT INTO result (set_id, pair_id, target, site, mask, multiplicity, failure)
                             SELECT ?, pair_id, target, site, mask, 1, NULL FROM v102''', (v102_id,))
        self.conn.commit()

    def v102_counts(self, target_name, mask):
        return { r[0] : r[1] for r in self.conn.execute('''SELECT v.site, count(v.rowid)
                                                           FROM v102 v
                                                           JOIN target t ON v.target=t.rowid
                                                           WHERE v.site != -1 AND t.name=? AND mask=?
                                                           GROUP BY site
                                                           ORDER BY site''', (target_name, mask,)) }

    def our_pairs_missing_from_v102(self):
        return [ ( r[0], r[1], r[2],
                   r[3], r[4], r[5] ) for r in self.conn.execute('''SELECT tr.name, r.site, r.mask, tv.name, v.site, v.mask
                                                                    FROM result r
                                                                    LEFT JOIN target tr ON r.target=tr.rowid
                                                                    LEFT JOIN v102 v ON r.pair_id = v.pair_id
                                                                    LEFT JOIN target tv ON v.target=tv.rowid
                                                                    WHERE r.site != -1 AND (v.site = -1 or v.rowid IS NULL)''') ]

    def v102_pairs_missing_from_ours(self):
        return [ ( r[0], r[1], r[2],
                   r[3], r[4], r[5] ) for r in self.conn.execute('''SELECT tr.name, r.site, r.mask, tv.name, v.site, v.mask
                                                                    FROM v102 v
                                                                    LEFT JOIN target tv ON v.target=tv.rowid
                                                                    LEFT JOIN result r ON r.pair_id = v.pair_id
                                                                    LEFT JOIN target tr ON r.target=tr.rowid
                                                                    WHERE v.site != -1 AND (r.site=-1 OR r.rowid IS NULL)''') ]

    def our_pairs_differing_from_v102(self):
        return [ ( r[0], r[1], r[2],
                   r[3], r[4], r[5] ) for r in self.conn.execute('''SELECT tr.name, r.site, r.mask, tv.name, v.site, v.mask
                                                                    FROM result r
                                                                    LEFT JOIN target tr ON r.target=tr.rowid
                                                                    JOIN v102 v ON r.pair_id = v.pair_id
                                                                    LEFT JOIN target tv ON v.target=tv.rowid
                                                                    WHERE r.site != -1 AND v.site != -1 AND v.site != r.site''') ]

    def our_pairs_nonmatching_v102(self):
        return [ ( r[0], r[1], r[2],
                   r[3], r[4], r[5] ) for r in self.conn.execute('''SELECT tr.name, r.site, r.mask, tv.name, v.site, v.mask
                                                                    FROM result r
                                                                    LEFT JOIN target tr ON r.target=tr.rowid
                                                                    JOIN v102 v ON r.pair_id = v.pair_id
                                                                    LEFT JOIN target tv ON v.target=tv.rowid
                                                                    WHERE r.site != v.site''') ]

    def v102_pairs_missing_from_ours_reasons(self):
        return [ ( r[0], r[1] ) for r in self.conn.execute('''SELECT r.failure, COUNT(v.pair_id)
                                                                    FROM v102 v
                                                                    LEFT JOIN result r ON r.pair_id = v.pair_id
                                                                    WHERE v.site != -1 AND (r.site=-1 OR r.rowid IS NULL)
                                                                    GROUP BY r.failure''') ]

    def our_pairs_nonmatching_v102_for_reason(self, reason):
        return self.conn.execute('''SELECT p.rowid, p.identifier, p.r1, p.r2, v.numeric_id, v.site, v.nomask_r1, v.nomask_r2
                                    FROM result r
                                    JOIN pair p ON p.rowid = r.pair_id
                                    JOIN v102 v ON r.pair_id = v.pair_id
                                    WHERE r.site != v.site AND r.failure=?''', (reason,))

    def diagram_info(self, rowid):
        return self.conn.execute('''SELECT p.rowid, p.identifier, p.r1, p.r2, v.numeric_id, v.site, v.nomask_r1, v.nomask_r2
                                    FROM result r
                                    JOIN pair p ON p.rowid = r.pair_id
                                    JOIN v102 v ON r.pair_id = v.pair_id
                                    WHERE r.pair_id=?''', (rowid,)).fetchone()

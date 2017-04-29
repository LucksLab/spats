
import os
import sqlite3

from parse import FastFastqParser

class PairDB(object):

    def __init__(self, database_path = None):
        self._dbpath = database_path or ":memory:"
        self.conn = sqlite3.connect(self._dbpath)

    def wipe(self):
        self.conn.execute("DROP TABLE IF EXISTS pair")
        self.conn.execute("DROP INDEX IF EXISTS r1_idx")
        self.conn.execute("DROP INDEX IF EXISTS r2_idx")
        self.conn.execute("DROP INDEX IF EXISTS pair_id_idx")
        self._wipe_v102()

    def _create(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS pair (r1 TEXT, r2 TEXT, identifier TEXT)")

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
        with FastFastqParser(r1_path, r2_path) as parser:
            while True:
                pairs, count = parser.read(16384)
                total += count
                if not pairs:
                    break
                conn.executemany("INSERT INTO pair (identifier, r1, r2) VALUES (?, ?, ?)", pairs)
        conn.commit()
        return total

    def index(self):
        self.conn.execute("CREATE INDEX IF NOT EXISTS r1_idx ON pair (r1)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS r2_idx ON pair (r2)")

    def _fetch_one(self, query, args = []):
        return self.conn.execute(query, args).fetchone()[0]

    def count(self):
        return self._fetch_one("SELECT COUNT(*) from pair")

    def unique_r1(self):
        self.index()
        return self._fetch_one("SELECT COUNT(distinct r1) from pair")

    def unique_r2(self):
        self.index()
        return self._fetch_one("SELECT COUNT(distinct r2) from pair")

    def unique_pairs(self):
        self.index()
        return self._fetch_one("SELECT COUNT(distinct r1||r2) from pair")

    def max_r1(self):
        self.index()
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r1 order by cnt desc limit 1")

    def max_r2(self):
        self.index()
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r2 order by cnt desc limit 1")

    def _batch_results(self, query, batch_size):
        batch = []
        count = 0
        for result in self.conn.execute(query):
            batch.append((int(result[0]), str(result[1]), str(result[2]), int(result[3])))
            count += 1
            if batch_size and count >= batch_size:
                cur = batch
                batch = []
                count = 0
                yield cur
        yield batch

    def unique_pairs_with_counts(self, batch_size = 0):
        return self._batch_results("SELECT count(rowid), r1, r2, rowid from pair group by (r1||r2)", batch_size)

    def all_pairs(self, batch_size = 0):
        print "Using all_pairs..."
        return self._batch_results("SELECT 1, r1, r2, rowid from pair", batch_size)

    def add_targets_table(self, targets_path):
        from parse import fasta_parse
        self.conn.execute("DROP TABLE IF EXISTS target")
        self.conn.execute("DROP INDEX IF EXISTS target_idx")
        self.conn.execute("CREATE TABLE target (name TEXT, seq TEXT)")
        self.conn.execute("CREATE INDEX target_idx ON target (name)")
        target_list = fasta_parse(targets_path)
        self.conn.executemany('INSERT INTO target (name, seq) VALUES (? , ?)', target_list)
        return { name : seq for name, seq in target_list }

    # results storing
    def prepare_results(self):
        self.conn.execute("DROP TABLE IF EXISTS result")
        self.conn.execute("CREATE TABLE IF NOT EXISTS result (pair_id INT, target INT, site INT, mask TEXT, multiplicity INT, failure TEXT)")
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS pair_result_idx ON result (pair_id)")

    # results a list of (rowid, target_name, site, mask, multiplicity)
    def add_results(self, results):
        # grab a new connection since this might be in a new process (due to multiprocessing)
        conn = sqlite3.connect(self._dbpath)
        cursor = conn.executemany('''INSERT INTO result (pair_id, target, site, mask, multiplicity, failure)
                                     VALUES (?, (SELECT rowid FROM target WHERE name=?), ?, ?, ?, ?)''', results)
        if cursor.rowcount != len(results):
            print results
            raise Exception("some results failed to add: {} / {}".format(cursor.rowcount, len(results)))
        conn.commit()

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


import os
import random
import sqlite3
import sys
import time

from parse import FastFastqParser

SHOW_SLOW_QUERIES = False

class _ConnWrapper(object):
    def __init__(self, conn):
        self._conn = conn
    def execute(self, query, args = []):
        start = time.time()
        res = self._conn.execute(query, args)
        delta = time.time() - start
        if delta > 1.0:
            print("VERY SLOW [{:.1f}s]: ".format(delta) + query)
        elif delta > 0.1:
            print("SLOW [{:.01f}s]: ".format(delta) + query)
        return res
    def executemany(self, query, multiargs):
        return self._conn.executemany(query, multiargs)
    def cursor(self):
        return self._conn.cursor()
    def commit(self):
        self._conn.commit()

class PairDB(object):

    def __init__(self, database_path = None):
        self._dbpath = database_path or ":memory:"
        self.conn = self._get_connection()
        self.show_progress_every = False

    def _get_connection(self):
        conn = sqlite3.connect(self._dbpath)
        return _ConnWrapper(conn) if SHOW_SLOW_QUERIES else conn

    def wipe(self):
        self.conn.execute("DROP TABLE IF EXISTS pair")
        self.conn.execute("DROP INDEX IF EXISTS r1_idx")
        self.conn.execute("DROP INDEX IF EXISTS r2_idx")
        self.conn.execute("DROP INDEX IF EXISTS pair_id_idx")
        self.conn.execute("DROP TABLE IF EXISTS unique_pair")
        self.conn.execute("DROP TABLE IF EXISTS tag")
        self.conn.execute("DROP TABLE IF EXISTS target")
        self.conn.execute("DROP TABLE IF EXISTS tag_count")
        self.conn.execute("DROP TABLE IF EXISTS result")
        self.conn.execute("DROP TABLE IF EXISTS result_set")
        self.conn.execute("DROP TABLE IF EXISTS result_tag")
        self._wipe_v102()

    def _create(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS pair (r1 TEXT, r2 TEXT, identifier TEXT)")

    def load_and_index(self, target_path, r1_path, r2_path):
        self.add_targets_table(target_path)
        self.parse(r1_path, r2_path)
        self._cache_unique_pairs()

    def _parse_and_sample(self, r1_path, r2_path, sample_size):
        conn = self.conn
        total = 0
        show_progress_every = self.show_progress_every
        next_report_count = show_progress_every
        with FastFastqParser(r1_path, r2_path) as parser:
            num_pairs_left = parser.appx_number_of_pairs()
            num_sampled = 0
            read_chunk_size = 65536
            print("Sampling {} out of ~{:.1f}M total pairs".format(sample_size, float(num_pairs_left)/1000000.0))
            while True:
                pairs, count = parser.read(read_chunk_size)
                if not pairs:
                    break
                if count < read_chunk_size:
                    num_to_take = sample_size - num_sampled
                else:
                    num_to_take = int((float(count) / float(max(1, num_pairs_left))) * float(sample_size - num_sampled))
                pairs = random.sample(pairs, min(num_to_take, count, sample_size - num_sampled))
                conn.executemany("INSERT INTO pair (identifier, r1, r2) VALUES (?, ?, ?)", pairs)
                num_sampled += len(pairs)
                num_pairs_left -= count
                total += count
                if show_progress_every:
                    next_report_count -= count
                    if next_report_count < 0:
                        next_report_count = show_progress_every
                        sys.stdout.write('.')
                        sys.stdout.flush()
        conn.commit()
        if show_progress_every:
            print(".")

        return total

    # we might be able to speed this up, but for now this seems fast enough (~300k pairs/s?)
    def parse(self, r1_path, r2_path, sample_size = 0):
        conn = self.conn
        conn.execute("DROP INDEX IF EXISTS r1_idx")
        conn.execute("DROP INDEX IF EXISTS r2_idx")
        self._create()

        if sample_size:
            return self._parse_and_sample(r1_path, r2_path, sample_size)

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
        if sample_size:
            samples = random.sample(xrange(total), min(sample_size, total))
            conn.execute("DELETE FROM pair WHERE (1+rowid) NOT IN (" + ",".join(map(str,samples)) + ")")
            conn.commit()
            conn.execute("VACUUM")
            total = sample_size
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

    def has_pairs(self):
        try:
            return (1 == self._fetch_one("SELECT 1 FROM pair LIMIT 1"))
        except:
            return False

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
        return self._batch_results("SELECT 1, r1, r2, rowid, rowid FROM pair WHERE rowid >= {} ORDER BY rowid LIMIT {}", batch_size)

    def _prep_targets(self):
        self.conn.execute("DROP TABLE IF EXISTS target")
        self.conn.execute("DROP INDEX IF EXISTS target_idx")
        self.conn.execute("CREATE TABLE target (name TEXT, seq TEXT)")
        self.conn.execute("CREATE INDEX target_idx ON target (name)")

    def add_targets_table(self, targets_path):
        self._prep_targets()
        from parse import fasta_parse
        target_list = fasta_parse(targets_path)
        self.conn.executemany('INSERT INTO target (name, seq) VALUES (? , ?)', target_list)
        self.conn.commit()
        return { name : seq for name, seq in target_list }

    def add_targets(self, targets):
        self._prep_targets()
        self.conn.executemany('INSERT INTO target (name, seq) VALUES (? , ?)', [ (t.name, t.seq) for t in targets.targets ])
        self.conn.commit()

    def targets(self):
        return [ (str(r[0]), str(r[1]), int(r[2])) for r in self.conn.execute("SELECT name, seq, rowid FROM target") ]

    # results storing
    def _prepare_results(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS result (set_id INT, pair_id INT, target INT, mask TEXT, site INT, end INT, mut INT, multiplicity INT, failure TEXT, tag_mask INT)")
        self.conn.execute("DROP INDEX IF EXISTS pair_result_idx")
        self.conn.execute("DROP INDEX IF EXISTS result_site_idx")
        self.conn.execute("CREATE TABLE IF NOT EXISTS result_set (set_id TEXT)")

    def add_result_set(self, set_name, resume_processing = False):
        self._prepare_results()
        rid = self.result_set_id_for_name(set_name)
        if rid is not None and not resume_processing:
            self.conn.execute("DELETE FROM result WHERE set_id=?", (rid,))
            self.conn.execute("DELETE FROM result_set WHERE rowid=?", (rid,))
            try:
                self.conn.execute("DELETE FROM result_tag WHERE set_id=?", (rid,))
            except:
                pass
            rid = None
        if rid is None:
            self.conn.execute("INSERT INTO result_set (set_id) VALUES (?)", (set_name,))
        self.conn.commit()
        return self.result_set_id_for_name(set_name)

    def result_set_id_for_name(self, set_name):
        self._prepare_results()
        return self._fetch_one("SELECT rowid FROM result_set WHERE set_id=?", (set_name,))

    def result_sets(self):
        self._prepare_results()
        return [ str(r[0]) for r in self.conn.execute("SELECT set_id FROM result_set") ]

    def add_results_with_tags(self, result_set_id, results):
        # grab a new connection since this might be in a new process (due to multiprocessing)
        conn = self._get_connection()
        stmt = '''INSERT INTO result (set_id, pair_id, target, mask, site, end, mut, multiplicity, failure, tag_mask)
                  VALUES ({}, ?, ?, ?, ?, ?, ?, ?, ?, ?)'''.format(result_set_id)
        rstmt_template = 'INSERT INTO result_tag (set_id, result_id, tag_id) VALUES ({}, {}, ?)'
        for res in results:
            cursor = conn.cursor()
            mask = sum([ 1 << t for t in res[-1] ])
            cursor.execute(stmt, [ mask if i == len(res) - 1 else res[i] for i in range(len(res)) ])
            rstmt = rstmt_template.format(result_set_id, cursor.lastrowid)
            cursor.executemany(rstmt, [ (t,) for t in res[-1] ])
        conn.commit()

    # results a list of (rowid, target_name, site, mask, multiplicity, failure, [optional: tags])
    def add_results(self, result_set_id, results):
        if len(results[0]) == 9:
            self.add_results_with_tags(result_set_id, results)
            return
        # grab a new connection since this might be in a new process (due to multiprocessing)
        #print " --> Thd in AR {}".format(self.worker_id)
        conn = sqlite3.connect(self._dbpath)
        stmt = '''INSERT INTO result (set_id, pair_id, target, mask, site, end, mut, multiplicity, failure)
                  VALUES ({}, ?, ?, ?, ?, ?, ?, ?, ?)'''.format(result_set_id)
        cursor = conn.executemany(stmt, results)
        if cursor.rowcount != len(results):
            raise Exception("some results failed to add: {} / {}".format(cursor.rowcount, len(results)))
        conn.commit()

    def index_results(self):
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS pair_result_idx ON result (set_id, pair_id)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS result_site_idx ON result (set_id, target, mask, end)")

    def num_results(self, result_set_name):
        return self._fetch_one("SELECT count(*) FROM result r JOIN result_set n ON r.set_id = n.rowid WHERE n.set_id = ?", (result_set_name,))

    def differing_results(self, result_set_name_1, result_set_name_2):
        id1 = self.result_set_id_for_name(result_set_name_1)
        id2 = self.result_set_id_for_name(result_set_name_2)
        return self.conn.execute('''SELECT p.rowid, p.r1, p.r2, s1.target, s1.end, s1.site, s1.mut, s1.failure, s2.target, s2.end, s2.site, s2.mut, s2.failure
                                    FROM result s1 JOIN result s2 ON s2.set_id=? AND s2.pair_id=s1.pair_id
                                    JOIN pair p ON p.rowid=s1.pair_id
                                    WHERE s1.set_id=? AND (s1.site != s2.site OR s1.end != s2.end OR s1.target != s2.target OR s1.mut != s2.mut) AND (s1.site != -1 OR s2.site != -1)''', (id2, id1))


    def result_sites(self, result_set_id, target_id):
        return self.conn.execute('''SELECT r.mask, r.end, r.site, SUM(r.multiplicity)
                                    FROM result r
                                    WHERE r.set_id = ? AND r.target = ? AND r.site != -1
                                    GROUP BY r.mask||r.end||r.site
                                    ORDER BY r.end, r.site ASC''', (result_set_id, target_id))

    # results analysis
    def setup_tags(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS tag (name TEXT)")
        self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS tag_name_idx ON tag (name)")
        self.conn.execute("CREATE TABLE IF NOT EXISTS result_tag (set_id INT, result_id INT, tag_id INT)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS result_tag_idx ON result_tag (set_id, tag_id)")

    def add_tags(self, tags):
        existing = self.tagmap()
        for t in tags:
            if -1 == existing.get(t, -1):
                self.conn.execute("INSERT INTO tag (name) VALUES (?)", (t,) )
        self.conn.commit()

    def tagmap(self):
        return { str(res[1]) : int(res[0]) for res in self.conn.execute("SELECT rowid, name FROM tag") }

    def count_tags(self, result_set_id):
        self.conn.execute("DROP TABLE IF EXISTS tag_count")
        self.conn.execute("CREATE TABLE tag_count (set_id INT, tag_id INT, count INT)")
        self.conn.execute('''INSERT INTO tag_count
                             SELECT ?, t.rowid, SUM(r.multiplicity)
                             FROM result_tag rt
                             JOIN tag t ON t.rowid=rt.tag_id
                             JOIN result r ON r.rowid=rt.result_id
                             GROUP BY t.rowid''', (result_set_id,) )
        self.conn.commit()

    def tag_counts(self, result_set_id, incl_tags = None, excl_tags = None):
        if not incl_tags and not excl_tags:
            query = "SELECT t.name, tc.count FROM tag_count tc JOIN tag t ON t.rowid = tc.tag_id WHERE set_id=?"
        else:
            query = self._tag_query(incl_tags, excl_tags)
            query = "SELECT t.name, SUM(mult) FROM ( " + query + " ) JOIN result_tag rt ON rrid=rt.result_id JOIN tag t ON t.rowid = rt.tag_id GROUP BY t.rowid"
        return { str(res[0]) : int(res[1]) for res in self.conn.execute(query, (result_set_id,)) }

    def _tag_query(self, incl_tags, excl_tags):
        tagmap = self.tagmap()
        tag_clause = ''
        def _mask(tags):
            return sum([ 1 << tagmap[t] for t in tags ])
        if incl_tags:
            mask = _mask(incl_tags)
            tag_clause += " AND (tag_mask & {} == {})".format(mask, mask)
        if excl_tags:
            mask = _mask(excl_tags)
            mask = sum([ 1 << tagmap[t] for t in excl_tags ])
            tag_clause += " AND (tag_mask & {} == 0)".format(mask)
        query = '''SELECT p.rowid, p.identifier, p.r1, p.r2, r.multiplicity as mult, r.rowid AS rrid, r.tag_mask
                   FROM result r
                   JOIN pair p ON p.rowid=r.pair_id
                   WHERE r.set_id=?''' + tag_clause + ' ORDER BY mult DESC'
        return query

    def results_matching(self, result_set_id, incl_tags = None, excl_tags = None, limit = 0):
        query = self._tag_query(incl_tags, excl_tags)
        if limit > 0:
            query += ' LIMIT {}'.format(limit)
        results =  self.conn.execute(query, (result_set_id,))
        return [ ( int(r[0]), str(r[1]), str(r[2]), str(r[3]), int(r[4]) ) for r in results ]

    def count_matches(self, result_set_id, incl_tags = None, excl_tags = None):
        query = self._tag_query(incl_tags, excl_tags)
        query = 'SELECT SUM(mult) FROM ( ' + query + ' )'
        return self._fetch_one(query, (result_set_id,))

    def results_matching_site(self, result_set_id, target_id, end, site, limit = 0):
        results =  self.conn.execute('''SELECT p.rowid, p.identifier, p.r1, p.r2, r.multiplicity
                                        FROM result r
                                        JOIN pair p ON p.rowid=r.pair_id
                                        WHERE r.set_id=? AND r.target=? AND r.site=? AND r.end=?
                                        ORDER BY r.multiplicity DESC ''' +
                                     'LIMIT {}'.format(limit) if limit > 0 else '',
                                     (result_set_id, target_id, site, end))
        return [ ( int(r[0]), str(r[1]), str(r[2]), str(r[3]), int(r[4]) ) for r in results ]

    def counter_data_for_results(self, result_set_id):
        results =  self.conn.execute('SELECT target, mask, site, end, multiplicity FROM result WHERE set_id=? AND site != -1 AND end != -1', (result_set_id,))
        return [ ( int(r[0]), str(r[1]), int(r[2]), int(r[3]), int(r[4]) ) for r in results ]


    # counters storage
    def setup_counters(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS counter (run_key TEXT, dict_index INT, count_key TEXT, count INT)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS counter_key_idx ON counter (run_key)")

    def has_counters(self):
        self.setup_counters()
        return (1 == self._fetch_one("SELECT 1 FROM counter LIMIT 1"))

    def store_counters(self, run_key, counters):
        self.setup_counters()
        self.conn.execute("DELETE FROM counter WHERE run_key=?", (run_key,))
        count_data = counters.count_data()
        row_data = [ (run_key, i, key, count_data[i][key]) for i in range(2) for key in count_data[i].keys() ]
        self.conn.executemany("INSERT INTO counter VALUES (?, ?, ?, ?)", row_data)
        self.conn.commit()

    def load_counters(self, run_key, counters):
        count_data = ( {}, {} )
        results =  self.conn.execute("SELECT dict_index, count_key, count FROM counter WHERE run_key=?", (run_key,))
        for r in results:
            count_data[int(r[0])][str(r[1])] = int(r[2])
        counters.reset()
        counters.update_with_count_data(count_data)
        return counters

    def setup_run(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS run_data (param_key TEXT, param_val TEXT)")

    def store_run(self, run):
        self.setup_run()
        self.conn.execute("DELETE FROM run_data")
        for k, v in run.config_dict().iteritems():
            self.conn.execute("INSERT INTO run_data VALUES (?, ?)", (k, str(v)))
        self.conn.commit()

    def load_run(self, run):
        config_dict = {}
        for r in self.conn.execute("SELECT param_key, param_val FROM run_data"):
            config_dict[str(r[0])] = str(r[1])
        run.load_from_config(config_dict)


    #v102 delta analysis
    def _create_v102(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS v102 (pair_id INT, numeric_id INT, target INT, site INT, end INT, nomask_r1 TEXT, nomask_r2 TEXT, mask TEXT)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS pair_id_idx ON pair (identifier)")

    def _wipe_v102(self):
        self.conn.execute("DROP TABLE IF EXISTS v102")
        self.conn.execute("DROP INDEX IF EXISTS v102_numeric_idx")
        self.conn.execute("DROP INDEX IF EXISTS pair_id_idx")
        self.conn.execute("DROP INDEX IF EXISTS v102_site_idx")

    def add_v102_comparison(self, targets_path, output_path, cotrans = False):
        import sys
        from parse import reactivities_parse, SamParser
        print("Creating index...")
        self.index()
        self._wipe_v102()
        self._create_v102()
        total = 0
        conn = self.conn

        # create the targets table
        targets = self.add_targets_table(targets_path)
        
        print("Parsing NOMASK*...")
        # parse NOMASK to get the ID<->numeric_id map
        with FastFastqParser(os.path.join(output_path, 'NOMASK_1.fq'), os.path.join(output_path, 'NOMASK_2.fq')) as parser:
            while True:
                pairs, count = parser.read_nomask(16384)
                total += count
                if not pairs:
                    break
                conn.executemany('''INSERT INTO v102
                                           (pair_id, numeric_id, target, site, end, nomask_r1, nomask_r2, mask)
                                    SELECT  pair.rowid,       ?,   NULL,   -1,  -1,         ?,         ?, NULL
                                    FROM pair WHERE pair.identifier=?''', pairs)
                sys.stdout.write('.')
                sys.stdout.flush()
        conn.commit()
        conn.execute("CREATE INDEX IF NOT EXISTS v102_numeric_idx ON v102 (numeric_id)")
        check = self._fetch_one("SELECT COUNT(*) from v102")
        if check != total:
            raise Exception("some v102 pairs did not match? {} != {}".format(check, total))
        print("\nParsed {} records from NOMASK_{{1,2}}.fq".format(total))

        # now, parse sam to get the numeric_id->site,end map
        total = 0
        mask_totals = {}
        masks = ('RRRY', 'YYYR')
        for mask in masks:
            mask_total = 0
            with SamParser(os.path.join(output_path, mask + '.sam'), targets) as parser:
                insert = 'UPDATE v102 SET target=' + ('1' if cotrans else '(SELECT rowid FROM target WHERE name=?)') + ', site=?, end=?, mask=? WHERE numeric_id=?'
                while True:
                    pairs, count = parser.read(16384, mask, cotrans = cotrans)
                    if not pairs:
                        break
                    total += count
                    mask_total += count
                    conn.executemany(insert, pairs)
                    sys.stdout.write('.')
                    sys.stdout.flush()
            mask_totals[mask] = mask_total
        conn.commit()
        conn.execute("CREATE INDEX IF NOT EXISTS v102_site_idx ON v102 (site)")
        print("\nAdded {} records from *.sam:".format(total))
        for mask in masks:
            print("  {}.sam: {}".format(mask, mask_totals[mask]))

        # now, compare against reactivities.out to verify our counts
        if cotrans:
            all_counts = { "{}::{}::{}".format(r[0],r[1],r[2]) : r[3] for r in self.conn.execute('''SELECT v.mask, v.site, v.end, count(v.rowid)
                                                                                                    FROM v102 v
                                                                                                    GROUP BY v.end, v.site, v.mask''') }
        else:
            all_counts = { "{}::{}::{}".format(r[0],r[1],r[2]) : r[3] for r in self.conn.execute('''SELECT t.name, v.mask, v.site, count(v.rowid)
                                                                                                    FROM v102 v
                                                                                                    JOIN target t ON v.target=t.rowid
                                                                                                    GROUP BY t.name, v.site, v.mask''') }

        count = 0
        for entry in reactivities_parse(os.path.join(output_path, 'reactivities.out')):
            # (target, rt_start, site, nuc, treated_count, untreated_count, beta, theta, c)
            if cotrans:
                key = "{}::{}::{}".format('RRRY', int(entry[2]), int(entry[1]) - 19)
            else:
                key = "{}::{}::{}".format(entry[0], 'RRRY', int(entry[2]))
            if all_counts.get(key, 0) != int(entry[4]):
                raise Exception("Treated db does not match reactivities: key={}, db={}, r={}".format(key, all_counts.get(key, 0), int(entry[4])))
            if cotrans:
                key = "{}::{}::{}".format('YYYR', int(entry[2]), int(entry[1]) - 19)
            else:
                key = "{}::{}::{}".format(entry[0], 'YYYR', int(entry[2]))
            if all_counts.get(key, 0) != int(entry[5]):
                raise Exception("Untreated db does not match reactivities: key={}, db={}, r={}".format(key, all_counts.get(key, 0), int(entry[5])))
            count += 1
        if not count:
            raise Exception("no sites found in reactivities.out?")
        print("reactivities.out check pass ({} sites).".format(count))
        self.add_v102_to_results()

    def add_v102_to_results(self):
        v102_id = self.add_result_set("v102")
        self.conn.execute('''INSERT INTO result (set_id, pair_id, target, site, end, mask, multiplicity, failure)
                             SELECT ?, pair_id, target, site, end, mask, 1, NULL FROM v102''', (v102_id,))
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

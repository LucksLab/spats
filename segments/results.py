import sqlite3

from .util import LoggingClass, M, fs


class ResultsDB(LoggingClass):

    def __init__(self, key):
        LoggingClass.__init__(self, initLogging = True)
        self.key = key
        self.conn = None
        self._commitInterval = (1 << 18)
        self._written = 0

    def create(self):
        if self.conn:
            return

        dbFile = "{}_results.db".format(self.key)
        conn = sqlite3.connect(dbFile)
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS results")
        cur.execute("DROP TABLE IF EXISTS reads")
        cols = [
            ("sequenceId", "integer"),
            ("success", "integer"),
            ("barcode", "string"),
            ("roi", "string"),
            ("overlap", "integer"),
            ("errorCount", "integer"),
        ]
        sql = "CREATE TABLE results ({})".format(",".join([ f"{c[0]} {c[1]}" for c in cols ]))
        cur.execute(sql)
        cols = [
            ("sequenceId", "integer"),
            ("mutType", "integer"),
            ("site", "integer"),
            ("size", "integer"),
            ("substSize", "integer"),
        ]
        sql = "CREATE TABLE reads ({})".format(",".join([ f"{c[0]} {c[1]}" for c in cols ]))
        cur.execute(sql)

        cur.execute("CREATE INDEX barcode_idx ON results ( barcode )")
        cur.execute("CREATE INDEX res_seqid_idx ON results ( sequenceId )")
        cur.execute("CREATE INDEX reads_seqid_idx ON reads ( sequenceId )")

        conn.commit()
        self.conn = conn

    def write(self, results):
        cur = self.conn.cursor()
        for fr in results:
            seqid = fr.r1.sequenceId
            cur.execute('INSERT INTO results VALUES (?, ?, ?, ?, ?, ?)', (seqid,
                                                                          fr.c.success if fr.c else 0,
                                                                          fr.barcode,
                                                                          fr.c.roi if hasattr(fr.c, 'roi') else None,
                                                                          fr.f.overlap,
                                                                          fr.c.errorCount if fr.c else 0))
            if fr.c.errorCount and fr.c.mutations:
                cur.executemany('INSERT INTO reads VALUES (?, ?, ?, 1, 0)',
                                [ [ seqid, fs.mutation, e ] for e in fr.c.mutations])

            if fr.c.errorCount and fr.c.indels:
                cur.executemany('INSERT INTO reads VALUES (?, ?, ?, ?, ?)',
                                [ [ seqid, i.mutType, i.errorIndex, i.size, i.substSize ] for i in fr.c.indels])
            self._written += 1
            if self._written >= self._commitInterval:
                self.conn.commit()
                self._written = 0

    def finish(self):
        self.conn.commit()
        self.conn.close()
        self.conn = None
        self.info("Wrote {}_results.db".format(self.key))

    def save(self, results):
        self.create()
        self.write(results)
        self.finish()

    def barcodeCounts(self, limit = None):
        dbFile = "{}_results.db".format(self.key)
        conn = sqlite3.connect(dbFile)
        cur = conn.cursor()
        bcSql = "SELECT barcode, count(sequenceId) as cnt FROM results WHERE success=1 GROUP BY barcode ORDER BY cnt DESC {}".format(f'LIMIT {limit}' if limit else '')
        counts = []
        for bres in cur.execute(bcSql).fetchall():
            counts.append(bres[1])
        return counts

    def barcodeAnalysisInfo(self, largeMutSize = 4, limit = None):
        dbFile = "{}_results.db".format(self.key)
        conn = sqlite3.connect(dbFile)
        cur = conn.cursor()

        def equivalentRead(r1, r2):
            if len(r1.mutations) != len(r2.mutations):
                return False
            for m1, m2 in zip(r1.mutations, r2.mutations):
                if (m1.site != m2.site) or (m1.mutType != m2.mutType) or (m1.size != m2.size) or (m1.substSize != m2.substSize):
                    return False
            return True

        def hasLargeMut(rd):
            for m in rd.mutations:
                if m.size > largeMutSize or m.substSize > largeMutSize:
                    return True
            return False

        def incrementCount(count, array):
            if len(array) <= count:
                array += ([0] * (count + 1 - len(array)))
            array[count] += 1

        bcSql = "SELECT DISTINCT(barcode) FROM results WHERE success=1 {}".format(f'LIMIT {limit}' if limit else '')
        readSql = '''SELECT r.sequenceId, rd.mutType, rd.site, rd.size, rd.substSize, r.roi
                     FROM results r
                     LEFT JOIN reads rd ON rd.sequenceId=r.sequenceId
                     WHERE r.success AND r.barcode='{}'
                     ORDER BY rd.sequenceId, rd.site'''

        atInfo = {}
        multiplicityInfo = {}
        bcCount = 0
        self.info("  fetching barcodes..")
        for bres in cur.execute(bcSql).fetchall():
            bcCount += 1
            if 0 == (bcCount % 10000):
                self.info("  bcInfo {}..".format(bcCount))
            barcode = bres[0]
            barcodeReads = {} # sequenceID -> read
            rois = set()
            for rres in cur.execute(readSql.format(barcode)).fetchall():
                assert(rres[0] != None)
                rd = barcodeReads.setdefault(rres[0], M.NGSRead(sequenceId = rres[0]))
                rois.add(rres[5])
                rd.barcode = barcode
                if rres[1]:
                    assert(rres[1] in [fs.mutation, fs.insertion, fs.deletion, fs.substitution])
                    rd.mutations.append(M.NGSMutation(mutType = rres[1], site = rres[2], size = rres[3], substSize = rres[4]))
            multiplicity = len(barcodeReads)
            self.debug("  {}: {}".format(barcode, len(barcodeReads)), barcodeReads)
            for chopt in [ 'A', 'C', 'G', 'T', 'AT', 'AC', 'AG', 'CG', 'CT', 'GT' ]:
                chInfo = atInfo.setdefault(chopt, {})
                chCount = len([ch for ch in barcode if ch in chopt])
                chInfo.setdefault(chCount, []).append(multiplicity)
            bci = multiplicityInfo.setdefault(multiplicity, M.BarcodeInfo(multiplicity = multiplicity))
            bci.numBarcodes += 1

            readEquivalenceClasses = []
            perfect = muted = largeMuted = 0
            for rd in barcodeReads.values():
                if not rd.mutations:
                    perfect += 1
                else:
                    muted +=1
                    if hasLargeMut(rd):
                        largeMuted += 1

                foundClass = False
                for rec in readEquivalenceClasses:
                    if equivalentRead(rd, rec[0]):
                        rec.append(rd)
                        foundClass = True
                        break
                if not foundClass:
                    readEquivalenceClasses.append([rd])

            self.debug(barcodeReads, readEquivalenceClasses)
            assert(readEquivalenceClasses)
            if 1 == len(readEquivalenceClasses) and (1 == len(rois)):
                bci.agree += 1
            else:
                bci.disagree +=1

            bci.largeMutCount += largeMuted
            assert(len(readEquivalenceClasses) <= multiplicity)

            maj = [ x for x in readEquivalenceClasses if (len(x) << 1) > len(barcodeReads) ]
            if maj:
                assert(1 == len(maj))
                maj = maj[0][0]
                if not maj.mutations:
                    bci.perfectMajority += 1
                else:
                    bci.mutedMajority += 1
                if hasLargeMut(maj):
                    bci.largeMutMajority += 1
                else:
                    bci.largeMutDropped += largeMuted
            else:
                bci.noMajority += 1
                bci.largeMutDropped += largeMuted

            incrementCount(perfect, bci.perfectMatchCounts)
            incrementCount(muted, bci.mutedCounts)
            incrementCount(largeMuted, bci.largeMutedCounts)
            incrementCount(len(readEquivalenceClasses), bci.distinctReadCounts)
            incrementCount(max([len(rec) for rec in readEquivalenceClasses]), bci.maxAgreeingReadCounts)
            incrementCount(max([0] + [len(rec) for rec in readEquivalenceClasses if rec[0].mutations]), bci.maxAgreeingMutedReadCounts)
            incrementCount(max([0] + [len(rec) for rec in readEquivalenceClasses if hasLargeMut(rec[0])]), bci.maxAgreeingLargeMutedReadCounts)

        return multiplicityInfo, atInfo

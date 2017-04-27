
import sqlite3

from parse import FastFastqParser

class PairDB(object):

    def __init__(self, database_path):
        self._dbpath = database_path
        self.conn = sqlite3.connect(self._dbpath)

    def wipe(self):
        self.conn.execute("DROP TABLE pair")

    def setup(self):
        self.conn.execute("CREATE TABLE IF NOT EXISTS pair (r1 TEXT, r2 TEXT, identifier TEXT)")

    # we might be able to speed this up by:
    #  - some sqlite tricks like: conn.execute("PRAGMA synchronous=OFF"), conn.execute("PRAGMA cache_size=16000"), etc.
    #  - parsing in one thread and inserting in another
    # but for now this seems fast enough (~300k pairs/s?)
    def parse(self, r1_path, r2_path):
        conn = self.conn
        conn.execute("DROP INDEX IF EXISTS r1_idx")
        conn.execute("DROP INDEX IF EXISTS r2_idx")
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
        conn = self.conn
        conn.execute("CREATE INDEX IF NOT EXISTS r1_idx ON pair (r1)")
        conn.execute("CREATE INDEX IF NOT EXISTS r2_idx ON pair (r2)")

    def _fetch_one(self, query, args = []):
        return self.conn.execute(query, args).fetchone()[0]

    def count(self):
        return self._fetch_one("SELECT COUNT(*) from pair")

    def unique_r1(self):
        return self._fetch_one("SELECT COUNT(distinct r1) from pair")

    def unique_r2(self):
        return self._fetch_one("SELECT COUNT(distinct r2) from pair")

    def max_r1(self):
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r1 order by cnt desc limit 1")

    def max_r2(self):
        return self._fetch_one("SELECT COUNT(rowid) as cnt from pair group by r2 order by cnt desc limit 1")

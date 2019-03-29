
import unittest

from spats_shape_seq.experiment import Experiment, FragmentDescriptor, SectionMatcher, Sequence, Target


class MockExperiment(Experiment):

    def __init__(self):
        Experiment.__init__(self, "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHI")
        self.linker = Sequence("01234567890123456789", Sequence.LINKER)
        self.dumbbell = Sequence("0123456789012345678901", Sequence.DUMBBELL)
        self.adapter_t = Sequence("abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdef", Sequence.ADAPTERT)
        self.adapter_b = Sequence("abcdefghijklmnopqrstuvwxyzabcdefgh", Sequence.ADAPTERB)


class TestMatcher(unittest.TestCase):

    def test_sections_of_lengths(self):
        t = Target("0123")
        m = SectionMatcher(t, True, SectionMatcher.MATCH_SUBSTRING)
        ss = m.sectionsOfLengths(1, 4)
        self.assertEqual(10, len(ss))


class TestFragmentDescriptor(unittest.TestCase):

    def setUp(self):
        self.experiment = MockExperiment()

    def tearDown(self):
        self.experiment = None

    def test_one_exact(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l, l)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        pf = fd.perfectFragments()
        self.assertEqual(1, len(pf))
        self.assertEqual(x.dumbbell.string, pf[0].string)

    def test_one_prefix(self):
        x = self.experiment
        fd = FragmentDescriptor(5, 20)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(16, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        self.assertTrue(db[:4] not in pfstrs)
        self.assertTrue(db[:5] in pfstrs)
        self.assertTrue(db[:6] in pfstrs)
        self.assertTrue(db[:8] in pfstrs)
        self.assertTrue(db[:19] in pfstrs)
        self.assertTrue(db[:20] in pfstrs)
        self.assertTrue(db[:21] not in pfstrs)

    def test_one_suffix(self):
        x = self.experiment
        fd = FragmentDescriptor(5, 20)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_SUFFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(16, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        self.assertTrue(db[-4:] not in pfstrs)
        self.assertTrue(db[-5:] in pfstrs)
        self.assertTrue(db[-6:] in pfstrs)
        self.assertTrue(db[-8:] in pfstrs)
        self.assertTrue(db[-19:] in pfstrs)
        self.assertTrue(db[-20:] in pfstrs)
        self.assertTrue(db[-21:] not in pfstrs)

    def test_small_substring(self):
        s = Sequence("01234", Sequence.DUMBBELL)
        fd = FragmentDescriptor(2, 4)
        fd.addMatcher(SectionMatcher(s, True, SectionMatcher.MATCH_SUBSTRING))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(9, len(pf))
        pfstrs = [ f.string for f in pf ]
        ss = s.string
        self.assertTrue(ss[:1] not in pfstrs)
        self.assertTrue(ss[:2] in pfstrs)
        self.assertTrue(ss[:3] in pfstrs)
        self.assertTrue(ss[:4] in pfstrs)
        self.assertTrue(ss[:5] not in pfstrs)
        self.assertTrue(ss[1:4] in pfstrs)
        self.assertTrue(ss[2:5] in pfstrs)
        self.assertTrue(ss[-2:] in pfstrs)

    def test_one_substring(self):
        x = self.experiment
        fd = FragmentDescriptor(5, 20)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_SUBSTRING))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(168, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        self.assertTrue(db[:4] not in pfstrs)
        self.assertTrue(db[:5] in pfstrs)
        self.assertTrue(db[:6] in pfstrs)
        self.assertTrue(db[:8] in pfstrs)
        self.assertTrue(db[:19] in pfstrs)
        self.assertTrue(db[:20] in pfstrs)
        self.assertTrue(db[1:6] in pfstrs)
        self.assertTrue(db[1:21] in pfstrs)
        self.assertTrue(db[3:8] in pfstrs)
        self.assertTrue(db[3:] in pfstrs)
        self.assertTrue(db[-5:] in pfstrs)
        self.assertTrue(db[-6:] in pfstrs)
        self.assertTrue(db[-8:] in pfstrs)
        self.assertTrue(db[-19:] in pfstrs)
        self.assertTrue(db[-20:] in pfstrs)

    def test_exact_prefix(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l, l + 10)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(11, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        ad = x.adapter_t.string
        self.assertTrue(db in pfstrs)
        self.assertTrue(db + ad[:1] in pfstrs)
        self.assertTrue(db + ad[1:3] not in pfstrs)
        self.assertTrue(db + ad[:4] in pfstrs)
        self.assertTrue(db + ad[:10] in pfstrs)
        self.assertTrue(db + ad[:11] not in pfstrs)

    def test_exact_prefix_req(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l, l + 10)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.adapter_t, True, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(10, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        ad = x.adapter_t.string
        self.assertTrue(db not in pfstrs)
        self.assertTrue(db + ad[:1] in pfstrs)
        self.assertTrue(db + ad[1:3] not in pfstrs)
        self.assertTrue(db + ad[:4] in pfstrs)
        self.assertTrue(db + ad[:10] in pfstrs)
        self.assertTrue(db + ad[:11] not in pfstrs)

    def test_exact_suffix(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l, l + 10)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.target, False, SectionMatcher.MATCH_SUFFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(11, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        t = x.target.string
        self.assertTrue(db in pfstrs)
        self.assertTrue(db + t[-1:] in pfstrs)
        self.assertTrue(db + t[-3:-1] not in pfstrs)
        self.assertTrue(db + t[-4:] in pfstrs)
        self.assertTrue(db + t[-10:] in pfstrs)
        self.assertTrue(db + t[-11:] not in pfstrs)

    def test_exact_suffix_req(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l, l + 10)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUFFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(10, len(pf))
        pfstrs = [ f.string for f in pf ]
        db = x.dumbbell.string
        t = x.target.string
        self.assertTrue(db not in pfstrs)
        self.assertTrue(db + t[-1:] in pfstrs)
        self.assertTrue(db + t[-3:-1] not in pfstrs)
        self.assertTrue(db + t[-4:] in pfstrs)
        self.assertTrue(db + t[-10:] in pfstrs)
        self.assertTrue(db + t[-11:] not in pfstrs)

    def test_short_suffix_prefix(self):
        x = self.experiment
        fd = FragmentDescriptor(2, 4)
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUFFIX))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(9, len(pf))
        pfstrs = [ f.string for f in pf ]
        t = x.target.string
        ad = x.adapter_t.string
        self.assertTrue(t[-2:] in pfstrs)
        self.assertTrue(t[-1:] + ad[:1] in pfstrs)
        self.assertTrue(t[-1:] + ad[:2] in pfstrs)
        self.assertTrue(ad[:2] not in pfstrs)
        self.assertTrue(t[-4:] in pfstrs)
        self.assertTrue(t[-3:] + ad[:1] in pfstrs)

    # i.e., single-length experiment
    def test_suffix_prefix(self):
        x = self.experiment
        fd = FragmentDescriptor(20, 22)
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUFFIX))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(63, len(pf))

    # i.e., single-length dumbbell experiment
    def test_exact_suffix_prefix(self):
        x = self.experiment
        l = x.dumbbell.len
        fd = FragmentDescriptor(l + 10, l + 12)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUFFIX))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(33, len(pf))
        db = x.dumbbell.string
        t = x.target.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        self.assertTrue(db + t[-13:] not in pfstrs)
        self.assertTrue(db + t[-12:] in pfstrs)
        self.assertTrue(db + t[-3:] + ad[:8] in pfstrs)
        self.assertTrue(db + t[-1:] + ad[:11] in pfstrs)
        self.assertTrue(db + ad[:11] not in pfstrs)

    def test_substr_exact_prefix_short(self):
        x = self.experiment
        x.target = Target("ABCDEF")
        x.linker = Sequence("0123", Sequence.LINKER)
        ll = x.linker.len
        fd = FragmentDescriptor(ll + 1, ll + x.target.len)
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUBSTRING))
        fd.addMatcher(SectionMatcher(x.linker, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        pf = fd.perfectFragments()
        print("\n".join(map(str,pf)))
        self.assertEqual(91, len(pf))
        t = x.target.string
        l = x.linker.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        self.assertTrue(t + l in pfstrs)
        self.assertTrue(t[2:4] + l + ad[:1] in pfstrs)
        self.assertTrue(t[2:3] + l + ad[:1] in pfstrs)

    # i.e., cotrans experiment
    def test_substr_exact_prefix(self):
        x = self.experiment
        ll = x.linker.len
        fd = FragmentDescriptor(ll + 12, ll + x.target.len)
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUBSTRING))
        fd.addMatcher(SectionMatcher(x.linker, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX, fillToMatchFragmentLength = 36))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(4126, len(pf))
        t = x.target.string
        l = x.linker.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        self.assertTrue(t + l in pfstrs)
        self.assertTrue(t[1:] + l in pfstrs)
        self.assertTrue(t[:-1] + l in pfstrs)
        self.assertTrue(t[8:11] + l + ad[:36-len(l)-3] in pfstrs)

    # i.e., cotrans+dumbbell experiment
    def test_exact_substr_exact_prefix(self):
        x = self.experiment
        dl = x.dumbbell.len
        ll = x.linker.len
        fd = FragmentDescriptor(dl + ll + 12, dl + ll + x.target.len)
        fd.addMatcher(SectionMatcher(x.dumbbell, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.target, True, SectionMatcher.MATCH_SUBSTRING))
        fd.addMatcher(SectionMatcher(x.linker, True, SectionMatcher.MATCH_EXACT))
        fd.addMatcher(SectionMatcher(x.adapter_t, False, SectionMatcher.MATCH_PREFIX, fillToMatchFragmentLength = 64))
        pf = fd.perfectFragments()
        #print("\n".join(map(str,pf)))
        self.assertEqual(4543, len(pf))
        d = x.dumbbell.string
        t = x.target.string
        l = x.linker.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        self.assertTrue(d + t + l in pfstrs)
        self.assertTrue(d + t[1:] + l in pfstrs)
        self.assertTrue(d + t[:-1] + l in pfstrs)
        self.assertTrue(d + t[8:11] + l + ad[:64-len(l)-len(d)-3] in pfstrs)


class TestExperiment(unittest.TestCase):

    def setUp(self):
        self.experiment = MockExperiment()

    def tearDown(self):
        self.experiment = None

    def test_single(self):
        x = self.experiment
        pf = x.descriptor().perfectFragments()
        self.assertEqual(91, len(pf))
        #print("\n".join(map(str,pf)))
        t = x.target.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        r1len = x.r1Length
        r2len = x.r2Length
        self.assertTrue(t in pfstrs)
        self.assertTrue(t[1:] in pfstrs)
        self.assertTrue(t[-r2len:] in pfstrs)
        self.assertTrue(t[-r2len:] + ad[:1] not in pfstrs)
        self.assertTrue(t[-r1len:] in pfstrs)
        self.assertTrue(t[-r1len:] + ad[:r2len-r1len] in pfstrs)
        self.assertTrue(t[-r1len + 1:] not in pfstrs)
        self.assertTrue(t[-r1len + 1:] + ad[:r2len-r1len+1] in pfstrs)
        self.assertTrue(t[-5:] + ad[:r2len - 5] in pfstrs)
        self.assertTrue(t[:-1] not in pfstrs)

    def test_single_dumbbell(self):
        x = self.experiment
        x.useDumbbell = True
        pf = x.descriptor().perfectFragments()
        self.assertEqual(91, len(pf))
        #print("\n".join(map(str,pf)))
        d = x.dumbbell.string
        t = x.target.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        r1len = x.r1Length
        r2len = x.r2Length
        dlen = x.dumbbell.len
        self.assertTrue(t not in pfstrs)
        self.assertTrue(d + t in pfstrs)
        self.assertTrue(d + t[1:] in pfstrs)
        self.assertTrue(d + t[-r2len:] in pfstrs)
        self.assertTrue(d + t[-r2len:] + ad[:1] not in pfstrs)
        self.assertTrue(d + t[-r1len:] in pfstrs)
        self.assertTrue(d + t[-(r1len-dlen):] + ad[:r2len-r1len] in pfstrs)
        self.assertTrue(d + t[-(r1len-dlen) + 1:] not in pfstrs)
        self.assertTrue(d + t[-(r1len-dlen)+1:] + ad[:r2len-r1len+1] in pfstrs)
        self.assertTrue(d + t[-5:] + ad[:(r2len - dlen) - 5] in pfstrs)
        self.assertTrue(d + t[:-1] not in pfstrs)

    def test_cotrans(self):
        x = self.experiment
        x.cotrans = True
        pf = x.descriptor().perfectFragments()
        self.assertEqual(4126, len(pf))
        #print("\n".join(map(str,pf)))
        l = x.linker.string
        t = x.target.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        r1len = x.r1Length
        r2len = x.r2Length
        llen = x.linker.len
        self.assertTrue(t not in pfstrs)
        self.assertTrue(t + l in pfstrs)
        self.assertTrue(t[1:] + l in pfstrs)
        self.assertTrue(t[7:86] + l in pfstrs)
        self.assertTrue(t[7:12] + l + ad[:r2len-llen-5] in pfstrs)
        self.assertTrue(t[:-1] + l in pfstrs)

    def test_cotrans_dumbbell(self):
        x = self.experiment
        x.cotrans = True
        x.useDumbbell = True
        pf = x.descriptor().perfectFragments()
        self.assertEqual(3828, len(pf))
        print("\n".join(map(str,pf)))
        l = x.linker.string
        d = x.dumbbell.string
        t = x.target.string
        ad = x.adapter_t.string
        pfstrs = [ f.string for f in pf ]
        r1len = x.r1Length
        r2len = x.r2Length
        llen = x.linker.len
        dlen = x.linker.len
        self.assertTrue(t not in pfstrs)
        self.assertTrue(t + l not in pfstrs)
        self.assertTrue(d + t not in pfstrs)
        self.assertTrue(d + l not in pfstrs)
        self.assertTrue(d + t + l in pfstrs)
        self.assertTrue(d + t[1:] + l in pfstrs)
        self.assertTrue(d + t[7:86] + l in pfstrs)
        self.assertTrue(d + t[7:12] + l in pfstrs)
        self.assertTrue(d + t[:-1] + l in pfstrs)


import cjb.uif
from cjb.uif.layout import Size, Grid, Rect
from cjb.uif.views import Label
from spats_shape_seq.util import reverse_complement
from viz.scenes.base import BaseScene
from viz.layout import buttonSize


class Nuc(object):

    def __init__(self, char, context):
        self.char = char
        self.context = context

    @property
    def displayName(self):
        return self.char

TAG_COLORS = [
    [ 0.7, 0.7, 0.7 ],
    [ 0.9, 0.6, 0.6 ],
    [ 0.6, 0.6, 0.9 ],
    [ 0.9, 0.4, 0.1 ],
    [ 0.2, 0.7, 0.3 ],
]

class RawPairScene(BaseScene):

    def __init__(self, ui, pair, expanded = False):
        self.pair = pair
        self.parts = {}
        self.labels = {}
        self.nucSize = Size(12, 18)
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addNucView(self, nuc, bg):
        v = cjb.uif.views.Button(obj = nuc)
        v.sideSpacing = 0
        v.bg = bg
        self.addView(v)
        return v

    def build(self):

        BaseScene.build(self)

        processor = self.ui.processor
        seqs = {
            "adapter_t_rc" : reverse_complement(processor._run.adapter_t),
            "adapter_b" : processor._run.adapter_b,
            "RRRY" : "RRRY",
            "YYYR" : "YYYR",
        }
        for target in processor._tag_targets.targets:
            seqs[target.name] = target.seq
            seqs[target.name + "_rc"] = reverse_complement(target.seq)
        self.tag_seqs = seqs

        tag_colors = {}
        colors = self.ui.colors
        for col in colors._colors.keys():
            tag_colors[col] = colors.color(col)
        color_idx = 1
        for part_name in ( "r1", "r2" ):
            part = getattr(self.pair, part_name)
            seq = part.original_seq
            self.parts[part_name] = map(lambda x: self.addNucView(x, TAG_COLORS[0]), [ Nuc(seq[idx], (part_name, idx, None)) for idx in range(len(seq)) ])
            self.labels[part_name] = self.addView(Label(part_name.upper()))
            for tindex in range(len(part.tags)):
                tag = part.tags[tindex]
                tkey = tag[0].rstrip("_rc_")
                tcol = tag_colors.get(tkey)
                if not tcol:
                    tcol = TAG_COLORS[color_idx % len(TAG_COLORS)]
                    color_idx += 1
                    tag_colors[tkey] = tcol
                tseq = self.tag_seqs[tag[0]]
                match_index = tag[3]
                views = []
                for idx in range(max(0,match_index - 4), min(len(tseq), match_index + tag[2] + 4)):
                    v = self.addNucView(Nuc(tseq[idx], (part_name, idx, tag)), tcol)
                    if idx < match_index or idx >= match_index + tag[2]:
                        v.alpha = 0.5
                        #print " set a=0.5, bg: {}".format(v.bg)
                    views.append(v)
                self.parts[part_name + str(tindex)] = views
                label = Label(tag[0])
                label.bg = tcol
                self.labels[part_name + str(tindex)] = self.addView(label)


    def layout(self, view):
        BaseScene.layout(self, view)
        cols = 100
        rows = 40
        frame = view.frame.centeredSubrect(self.nucSize.w * cols, self.nucSize.h * rows)
        grid = Grid(frame = frame, itemSize = self.nucSize, columns = cols, rows = rows)
        row_idx = 4

        def labelFrame(owner_parts, tag_parts):
            base_frame = tag_parts[0].frame.copy()
            base_frame.origin.x = owner_parts[0].frame.origin.x - 130
            base_frame.size = Size(80, base_frame.size.height)
            return base_frame

        for part_name in ( "r1", "r2" ):
            part_start = 32
            part = getattr(self.pair, part_name)
            if part_name == "r2":
                row_idx = 20
            grid.setLocation(part_start, row_idx)
            grid.applyToViews(self.parts[part_name])
            self.labels[part_name].frame = labelFrame(self.parts[part_name], self.parts[part_name])

            for tindex in range(len(part.tags)):
                tag = part.tags[tindex]
                tseq = self.tag_seqs[tag[0]]
                match_index = tag[3]
                tstart = max(0, match_index - 4)
                col = part_start + tag[1] - (match_index - tstart)
                row_idx += 1
                #print part_name + tag[0] + ": {} , {}".format(col, row_idx)
                grid.setLocation(col, row_idx)
                grid.applyToViews(self.parts[part_name + str(tindex)])

                self.labels[part_name + str(tindex)].frame = labelFrame(self.parts[part_name], self.parts[part_name + str(tindex)])

        return view


class PairInTargetScene(BaseScene):

    def __init__(self, ui, pair, expanded = True):
        self.pair = pair
        self.expanded = expanded
        self.parts = {}
        self.labels = {}
        self.nucSize = Size(10, 16)
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addNucView(self, nuc, bg):
        v = cjb.uif.views.Button(obj = nuc)
        v.fontSize = 11
        v.sideSpacing = 0
        v.bg = bg
        self.addView(v)
        return v

    def build(self):

        BaseScene.build(self)

        processor = self.ui.processor
        seqs = {
            "adapter_t_rc" : reverse_complement(processor._run.adapter_t),
            "adapter_b" : processor._run.adapter_b,
            "RRRY" : "RRRY",
            "YYYR" : "YYYR",
        }
        for target in processor._tag_targets.targets:
            seqs[target.name] = target.seq
            seqs[target.name + "_rc"] = reverse_complement(target.seq)
        self.tag_seqs = seqs

        colors = self.ui.colors
        target = self.pair.target
        colors._colors[target.name] = colors.color("target")
        tcol = colors.color("target")
        nomatch_col = colors.color("grey")
        error_col = colors.color("error")

        skips = self._skips()

        tseq = target.seq
        tlen = target.n
        def should_skip(idx):
            for skip in skips:
                if idx >= skip[0] and idx < skip[0] + skip[1]:
                    return True
            return False
        self.parts[target.name] = [ None if should_skip(i) else self.addNucView(Nuc(tseq[i], (target.name, i, None)), tcol) for i in range(tlen)  ]
        self.labels[target.name] = self.addLabel(target.name, bg = tcol)

        for part_name in ( "r1", "r2" ):
            parts = []
            part = getattr(self.pair, part_name)
            seq = part.original_seq
            idx = 0

            for tag in sorted(part.tags, key = lambda t : t[1]):
                tkey = tag[0].rstrip("_rc_")
                while idx < tag[1] + tag[2]:
                    ntcol = nomatch_col if idx < tag[1] else colors.color(tkey)
                    if idx in part.match_errors or idx in part.adapter_errors:
                        ntcol = error_col
                    parts.append(self.addNucView(Nuc(seq[idx], (part_name, idx, None if idx < tag[1] else tag[0])), ntcol))
                    idx += 1

                if self.expanded:
                    if tag[0] == target.name + "_rc":
                        rc = reverse_complement(seq[tag[1]:tag[1] + tag[2]])
                        self.parts[part_name + tag[0]] = [ self.addNucView(Nuc(rc[j], (tag[0], j, None)), tcol) for j in range(0, tag[2]) ]
                        self.labels[part_name + tag[0]] = self.addLabel("R1_rc")
                    elif tkey.startswith("adapter"):
                        tagseq = self.tag_seqs[tag[0]]
                        aparts = []
                        for j in range(max(tag[3] - 4, 0), min(tag[2] + 4, len(tagseq))):
                            v = self.addNucView(Nuc(tagseq[j], (tag[0], j, None)), colors.color(tkey))
                            if j < tag[3] or j >= tag[3] + tag[2]:
                                v.alpha = 0.5
                            aparts.append(v)
                        self.parts[part_name + tag[0]] = aparts
                        self.labels[part_name + tag[0]] = self.addLabel(tag[0], bg = colors.color(tkey))
                    elif part_name == "r1" and (tag[0] == 'YYYR' or tag[0] == 'RRRY'):
                        hcol = colors.color(tag[0])
                        self.parts[part_name + tag[0]] = [ self.addNucView(Nuc(tag[0][j], (tag[0], j, None)), hcol) for j in range(len(tag[0])) ]
                        self.labels[part_name + tag[0]] = self.addLabel(tag[0], bg = hcol)

            while idx < len(seq):
                parts.append(self.addNucView(Nuc(seq[idx], (part_name, idx, None)), nomatch_col))
                idx += 1

            self.parts[part_name] = parts
            self.labels[part_name] = self.addLabel(part_name.upper())


    def _skips(self):
        target = self.pair.target
        tlen = target.n
        r1matches = [ (tlen - tag[3] - tag[2], tag[2]) for tag in self.pair.r1.tags if tag[0].startswith(target.name) ] # tlen - b/c of revcomp
        r2matches = [ (tag[3], tag[2]) for tag in self.pair.r2.tags if tag[0].startswith(target.name) ]
        matches = sorted(r1matches + r2matches, key = lambda x : x[0])
        skips = []
        farthest = 0
        for idx in range(len(matches) + 1):
            m = matches[idx] if idx < len(matches) else None
            mp = matches[idx - 1] if idx > 0 else None
            if not mp:
                if m and m[0] > 0:
                    skips.append( (0, m[0]) )
            elif not m:
                if mp[0] + mp[1] < tlen:
                    skips.append( (mp[0] + mp[1], tlen - mp[0] - mp[1]) )
            else:
                right = max(farthest, mp[0] + mp[1])
                if right < m[0]:
                    skips.append( (right, m[0] - right) )
            if mp:
                farthest = max(farthest, mp[0] + mp[1])
        return [ (skip[0] + 4, skip[1] - 8) for skip in skips if skip[1] > 20 ]

    def layout(self, view):
        BaseScene.layout(self, view)
        cols = 100
        rows = 40
        frame = view.frame.centeredSubrect(self.nucSize.w * cols, self.nucSize.h * rows)
        grid = Grid(frame = frame, itemSize = self.nucSize, columns = cols, rows = rows)

        skips = self._skips()
        def skipped(idx):
            ret_idx = idx
            for skip in skips:
                if idx > skip[0]:
                    ret_idx = ret_idx - skip[1] + 1 # +1 for spaces representing skips
            return ret_idx

        target = self.pair.target
        row_idx = 4
        skipped_len = skipped(target.n)
        start_col = int((cols - skipped_len) / 2)

        def labelFrame(tag_parts):
            return Rect(grid.frame(start_col).origin.x - 140, tag_parts[0].frame.origin.y, 80, tag_parts[0].frame.size.height)

        grid.setLocation(start_col, row_idx)
        curskip_idx = 0
        nucs = self.parts[target.name]
        for idx in range(len(nucs)):
            nv = nucs[idx]
            if nv:
                nv.frame = grid.nextFrame()
            elif curskip_idx < len(skips) and idx == skips[curskip_idx][0]:
                nv = self.addNucView(Nuc("/", None), [1.0, 1.0, 1.0, 0.0 ])
                nv.frame = grid.nextFrame() # skip a frame
                curskip_idx += 1
        self.labels[target.name].frame = labelFrame(self.parts[target.name])

        row_idx += (2 if self.expanded else 1)
        r2tagmap = { tag[0] : tag for tag in self.pair.r2.tags }
        r2match = r2tagmap.get(target.name, ("", 0, 0, 0))
        grid.setLocation(start_col + skipped(r2match[3] - r2match[1]), row_idx)
        grid.applyToViews(self.parts["r2"])
        self.labels["r2"].frame = labelFrame(self.parts["r2"])

        nucs = self.parts.get("r2" + "adapter_t_rc")
        if nucs:
            row_idx += 1
            atag = r2tagmap["adapter_t_rc"]
            grid.setLocation(start_col + skipped(r2match[3] - r2match[1] + atag[1] - atag[3]), row_idx)
            grid.applyToViews(nucs)
            self.labels["r2" + "adapter_t_rc"].frame = labelFrame(nucs)

        r1tagmap = { tag[0] : tag for tag in self.pair.r1.tags }
        r1match = r1tagmap.get(target.name + "_rc", ("", 0, 0, 0))
        r1start = target.n - r1match[3] - r1match[2]
        nucs = self.parts.get("r1" + target.name + "_rc")
        if nucs:
            row_idx += 1
            grid.setLocation(start_col + skipped(r1start), row_idx)
            grid.applyToViews(nucs)
            self.labels["r1" + target.name + "_rc"].frame = labelFrame(nucs)

        same_row = (not self.expanded) and (r1start >= r2match[3] + r2match[2])
        row_idx += (0 if same_row else 1)
        grid.setLocation(start_col + skipped(r1start - r1match[1]), row_idx)
        grid.applyToViews(self.parts["r1"])
        if same_row:
            self.labels["r2"].text = "R2 / R1"
            self.labels["r1"].alpha = 0
        else:
            self.labels["r1"].frame = labelFrame(self.parts["r1"])

        nucs = self.parts.get("r1" + "adapter_b")
        if nucs:
            row_idx += 1
            atag = r1tagmap["adapter_b"]
            grid.setLocation(start_col + skipped(r1start - r1match[1] + atag[1] - atag[3]), row_idx)
            grid.applyToViews(nucs)
            self.labels["r1" + "adapter_b"].frame = labelFrame(nucs)

        for handle in [ "YYYR", "RRRY" ]:
            nucs = self.parts.get("r1" + handle)
            if nucs:
                row_idx += 1
                grid.setLocation(start_col + skipped(r1start - r1match[1]), row_idx)
                grid.applyToViews(nucs)
                self.labels["r1" + handle].frame = labelFrame(nucs)
            
        return view

    def handleKeyEvent(self, keyInfo):
        handler = { "x" : self.expand }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def expand(self, message = None):
        self.ui.setScene(PairInTargetScene(self.ui, self.pair, expanded = not self.expanded))


PairScene = PairInTargetScene

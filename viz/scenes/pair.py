
import cjb.uif
from cjb.uif.layout import Size, Grid
from cjb.uif.views import Label
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

class PairScene(BaseScene):

    def __init__(self, ui, pair):
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
        # TODO
        self.tag_seqs = { 
            "5s" : "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC",
            "5s_rc" : "GTCCTTGGTGCCCGAGTCAGATGCCTGGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACGGCCGCCAGGCATCC",
            "adapter_t_rc" : "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
            "adapter_b" : "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "RRRY" : "RRRY",
            "YYYR" : "YYYR",
        }

        BaseScene.build(self)
        tag_colors = {}
        color_idx = 1
        for part_name in ( "r1", "r2" ):
            part = getattr(self.pair, part_name)
            seq = part.original_seq
            self.parts[part_name] = map(lambda x: self.addNucView(x, TAG_COLORS[0]), [ Nuc(seq[idx], (part_name, idx, None)) for idx in range(len(seq)) ])
            self.labels[part_name] = self.addView(Label(part_name.upper()))
            for tag in part.tags:
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
                self.parts[part_name + tag[0]] = views
                label = Label(tag[0])
                label.bg = tcol
                self.labels[part_name + tag[0]] = self.addView(label)


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

            for tag in part.tags:
                tseq = self.tag_seqs[tag[0]]
                match_index = tag[3]
                tstart = max(0, match_index - 4)
                col = part_start + tag[1] - (match_index - tstart)
                row_idx += 1
                #print part_name + tag[0] + ": {} , {}".format(col, row_idx)
                grid.setLocation(col, row_idx)
                grid.applyToViews(self.parts[part_name + tag[0]])

                self.labels[part_name + tag[0]].frame = labelFrame(self.parts[part_name], self.parts[part_name + tag[0]])

        return view

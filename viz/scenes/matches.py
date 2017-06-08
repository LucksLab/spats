
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from viz.scenes.base import BaseScene
from viz.scenes.pair import PairScene, RawPairScene
from viz.layout import buttonSize

from spats_shape_seq.pair import Pair



class Nuc(object):

    def __init__(self, char, idx, context):
        self.char = char
        self.idx = idx
        self.context = context

    @property
    def displayName(self):
        return self.char


class MatchedPair(object):

    def __init__(self, r1, r2, multiplicity, rowid, identifier):
        self.r1 = r1
        self.r2 = r2
        self.r1_nucs = [ Nuc(r1[idx], idx, self) for idx in range(len(r1)) ]
        self.r2_nucs = [ Nuc(r2[idx], idx, self) for idx in range(len(r2)) ]
        self.multiplicity = multiplicity
        self.rowid = rowid
        self.identifier = identifier

    def display_id(self):
        return ":".join(self.identifier.split(':')[-2:])

    #@property
    #def displayName(self):
    #    return self.r1 + "   " + self.r2


class Matches(BaseScene):

    def __init__(self, ui, tag):
        self.include_tags = [tag]
        self.exclude_tags = []
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addMatchView(self, matched_pair):
        v = cjb.uif.views.View(obj = matched_pair)
        v.name_label = v.addSubview(cjb.uif.views.Label(str(matched_pair.display_id()), fontSize = 11))
        v.mult_label = v.addSubview(cjb.uif.views.Label(str(matched_pair.multiplicity), fontSize = 11))
        v.mult_label.alignment = "right"
        v.bg = [ 0.8, 0.8, 0.8 ]

        tagged_pair = self.processed_pair(matched_pair)

        bg = [0.7, 0.7, 0.7]

        # TODO: share
        tcol = [ 1.0, 0.85, 0.7 ]
        target = tagged_pair.target
        tag_colors = {
            "5S" : tcol,  ## TODO
            "adapter_t" : [ 1.0, 0.5, 0.5 ],
            "adapter_b" : [ 0.5, 0.5, 1.0 ],
            "RRRY" : [ 0.2, 1.0, 0.1 ],
            "YYYR" : [ 0.2, 1.0, 0.1 ],
            "RYYY" : [ 0.2, 1.0, 0.1 ],
            "YRRR" : [ 0.2, 1.0, 0.1 ],
        }
        nomatch_col = [ 0.7, 0.7, 0.7 ]
        error_col = [ 0.9, 0.1, 0.2 ]

        v.r1_nucs = []
        for nuc in matched_pair.r1_nucs:
            tag = None
            bg = nomatch_col
            for tag in tagged_pair.r1.tags:
                if nuc.idx >= tag[1] and nuc.idx < tag[1] + tag[2]:
                    bg = tag_colors[tag[0].rstrip("_rc")]
                    break
            if nuc.idx in tagged_pair.r1.match_errors or nuc.idx in tagged_pair.r1.adapter_errors:
                bg = error_col
            nv = cjb.uif.views.Button(obj = nuc)
            nv.sideSpacing = 0
            nv.bg = bg
            v.addSubview(nv)
            self.addView(nv)
            v.r1_nucs.append(nv)

        v.r2_nucs = []
        for nuc in matched_pair.r2_nucs:
            tag = None
            bg = nomatch_col
            for tag in tagged_pair.r2.tags:
                if nuc.idx >= tag[1] and nuc.idx < tag[1] + tag[2]:
                    bg = tag_colors[tag[0].rstrip("_rc")]
                    break
            if nuc.idx in tagged_pair.r2.match_errors or nuc.idx in tagged_pair.r2.adapter_errors:
                bg = error_col
            nv = cjb.uif.views.Button(obj = nuc)
            nv.sideSpacing = 0
            nv.bg = bg
            v.addSubview(nv)
            self.addView(nv)
            v.r2_nucs.append(nv)

        v.target = matched_pair
        #v.addSubview(cjb.uif.views.Button(obj = matched_pair))
        #v = cjb.uif.views.Button(obj = matched_pair)
        #v.fontSize = 11
        self.addView(v)
        return v

    def build(self):
        BaseScene.build(self)
        pairs = self.ui.db.results_matching(self.ui.result_set_id, self.include_tags, self.exclude_tags, limit = 50)
        matches = [ MatchedPair(p[2], p[3], p[4], p[0], p[1]) for p in pairs ]
        self.matchViews = [ self.addMatchView(m) for m in matches ]

    def layoutMatch(self, view):
        grid = Grid(frame = view.frame.bounds(), itemSize = Size(10, 16), columns = 100, rows = 1)
        grid.setLocation(12, 0)
        grid.applyToViews(view.r1_nucs)
        grid.setLocation(50, 0)
        grid.applyToViews(view.r2_nucs)
        f = grid.frame(0)
        f.update(origin = f.origin, w = f.size.width * 10, h = f.size.height)
        view.name_label.frame = f
        f = grid.frame(89)
        f.update(origin = f.origin, w = f.size.width * 10, h = f.size.height)
        view.mult_label.frame = f

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 1000, h = view.frame.size.h - 100)
        self.scroller = layoutInScroller(self.matchViews, cur, Size(1000, 14), 2, self.scroller)
        for v in self.matchViews:
            self.layoutMatch(v)
        return view

    def processed_pair(self, matched_pair):
        pair = Pair()
        pair.set_from_data(matched_pair.identifier, matched_pair.r1, matched_pair.r2, matched_pair.multiplicity)
        self.ui.processor.process_pair_detail(pair)
        return pair
        
    def handleViewMessage(self, scene, obj, message):
        if obj and (isinstance(obj, MatchedPair) or isinstance(obj, Nuc)):
            mp = obj if isinstance(obj, MatchedPair) else obj.context
            pair = self.processed_pair(mp)
            if pair.has_site:
                self.ui.pushScene(PairScene(self.ui, pair, expanded = True))
            else:
                self.ui.pushScene(RawPairScene(self.ui, pair, expanded = True))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)


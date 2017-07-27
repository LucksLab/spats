
import copy

from cjb.uif.layout import Size
from cjb.uif.views import Button, Scroller


buttonSize = Size(120, 40)
itemSize = Size(200, 60)


def basic(view, scene):
    cur = view.frame.centeredSubrect(w = itemSize.w, h = view.frame.size.h)
    container = Scroller()
    container.frame = cur
    view.addSubview(container)
    cur = cur.bounds()
    for sv in copy.copy(view.subviews):
        if sv.obj:
            sv.frame = cur.topSubrect(itemSize.h, margin = 20)
            container.addSubview(sv)
            cur = cur.leftover
        elif isinstance(sv, Button) and sv.key != 'home':
            sv.frame = view.frame.topLeftSubrect(size = buttonSize, margin = 20)
    return view

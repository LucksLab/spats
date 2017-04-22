
def sp(n):
    return " " * n

class Diagram(object):

    def __init__(self, target, pair):
        self.target = target
        self.pair = pair
        self.prefix_len = 8
        self.bars = []

    def _make_prefix(self, label):
        bit = label[:min(5, len(label))]
        return bit + sp(self.prefix_len - len(bit))

    def _make_target_lines(self):
        d = self._make_prefix("site")
        for i in range(0, self.target.n, 10):
            marker = str(i)
            if i > 0:
                marker = marker.replace('0','o')
            d += marker
            d += sp(10 - len(marker))
        self._add_line(d)

        d = self._make_prefix(self.target.name)
        seq = self.target.seq.lower()
        for part in [ self.pair.r1, self.pair.r2 ]:
            if part.match_len:
                start = part.match_index
                end = start + part.match_len
                seq = seq[:start] + seq[start:end].upper() + seq[end:]
        d += seq
        self._add_line(d)

    def _make_part(self, part):
        is_R1 = (part == self.pair.r1)
        match_index = part.match_index if part.matched else (((self.target.n - part.length) >> 1) + (40 if is_R1 else -40))
        hdr = sp(match_index - (5 if is_R1 else 0))
        if is_R1:
            if self.pair.mask:
                hdr += self.pair.mask.chars
            else:
                hdr += "????"
            hdr += " "
        spacer = (part.length >> 1) - 1 - (2 if is_R1 else 0)
        hdr += sp(spacer)
        hdr += "R1" if is_R1 else "R2"
        hdr += sp(self.target.n - len(hdr))
        self._add_line(self._make_prefix("") + hdr)

        spaces = match_index
        if not part.matched:
            spaces -= (5 if is_R1 else 0)
        elif is_R1:
            spaces -= (part.length - part.match_len + 1) # subtract len b/c match is in terms of revcomp, +1 for .
        else:
            spaces -= part.match_start
        d = sp(spaces)
        d += (part.seq[:4] + ("." if is_R1 else "") + part.seq[4:])
        d += sp(self.target.n - len(d))
        d = self._make_prefix("R1" if is_R1 else "R2") + d
        self._add_line(d)
        return [ self.prefix_len + part.match_index, self.prefix_len + part.match_index + part.match_len - 1 ] if part.matched else []

    def _make_r1_rev(self):
        r1 = self.pair.r1
        rev = r1.reverse_complement 
        d = "(revcomp)"
        d += sp(r1.match_index + self.prefix_len - len(d))
        d += rev[r1.match_start:r1.match_start+r1.match_len]
        d += sp(self.target.n - len(d))
        self._add_line(d)

    def _make_result(self, part):
        d = sp(part.match_index)
        l = "l={}".format(part.match_index)
        r = "r={}".format(part.match_index + part.match_len)
        leftover = part.match_len - len(l) - len(r) - 4
        if leftover < 0:
            raise Exception("not enough room: {}".format(part.match_len))
        halfish = (leftover >> 1)
        bit = "^{}{}, {}{}^".format("-"*halfish,l,r,"-"*(leftover - halfish))
        d += bit
        d += sp(self.target.n - len(d))
        self._add_line(sp(self.prefix_len) + d)

    def _add_line(self, line):
        for bar_list in self.bars:
            for bar in bar_list:
                if len(line) > bar and line[bar] == ' ':
                    line = line[:bar] + '|' + line[(bar+1):]
        if line.startswith("@") or line.startswith("\\"):
            pass
        elif line.startswith(" "):
            line = "|" + line
        else:
            line = "+" + line
        self.lines.append(line)

    def _make_summary(self):
        result = "\===>  "
        if self.pair.has_site:
            result += "++SITE {}:{}".format(self.pair.mask.chars, self.pair.site)
        elif not self.pair.mask:
            result += "mask failure"
        elif not self.pair.r1.matched:
            result += "R1 match not found"
        elif not self.pair.r2.matched:
            result += "R2 match not found"
        else:
            result += "FAIL"
        self.lines.append(result)

    def make(self):
        self.lines = [ ]
        base_len = self.target.n + self.prefix_len

        self._add_line("@" + self.pair.identifier)

        r2_bars = self._make_part(self.pair.r2)
        self.bars.append(r2_bars)

        r1_bars = self._make_part(self.pair.r1)
        self.bars.append(r1_bars)

        if self.pair.mask:
            self._add_line(sp(base_len))
            self._make_r1_rev()

        self._add_line(sp(base_len))
        self._add_line(sp(base_len))

        self._make_target_lines()

        if self.pair.r2.matched:
            self._make_result(self.pair.r2)
            self.bars.remove(r2_bars)
        if self.pair.r1.matched:
            self._make_result(self.pair.r1)
        else:
            self._add_line(sp(base_len))

        self._make_summary()

        return "\n".join(self.lines)

def diagram(target, pair):
    return Diagram(target, pair).make()

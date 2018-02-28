
import string

from target import _Target
from util import reverse_complement


indeterminate_translator = string.maketrans("ABCDEFGHIJKLMNOPQRSTUVWXYZ"," ! !!! !!!!!!!!!!!! !!!!!!")

def sp(n, bit = " "):
    return bit * n

class Diagram(object):

    def __init__(self, pair, run):
        self.target = pair.target or _Target('???', '?' * 100, 0)
        self.pair = pair
        self.run = run
        self.prefix_len = 8
        self.bars = []
        self.max_len = 0

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
        match_index = part.match_index if part.matched else (((self.target.n - part.original_len) >> 1) + (40 if is_R1 else -40))
        hdr = sp(match_index - (5 if is_R1 else 0))
        if is_R1:
            if self.pair.mask:
                hdr += self.pair.mask.chars
            else:
                hdr += "????"
            hdr += " "
        spacer = (part.original_len >> 1) - 1 - (2 if is_R1 else 0)
        hdr += sp(spacer)
        hdr += "R1" if is_R1 else "R2"
        self._add_line(self._make_prefix("") + hdr)

        spaces = match_index - (5 if is_R1 else 0)
        d = sp(spaces)
        if is_R1:
            d += (part.original_seq[:4] + ".")
        d += part.subsequence
        if part._rtrim:
            trimmed = part.original_seq[-part._rtrim:]
            if is_R1 or len(trimmed) < 4:
                d += ("." + trimmed)
            else:
                d += ("." + trimmed[:4] + "." + trimmed[4:])
        d = self._make_prefix("R1" if is_R1 else "R2") + d
        self._add_line(d)

        if part.matched:
            return [ self.prefix_len + part.match_index, self.prefix_len + part.match_index + part.match_len - 1 ]
        elif self.pair.failure == "indeterminate sequence failure":
            d = sp(spaces)
            if is_R1:
                d += part.original_seq[:4].translate(indeterminate_translator) + "." + part.original_seq[4:].translate(indeterminate_translator)
            else:
                d += part.original_seq.translate(indeterminate_translator)
            self._add_line(self._make_prefix("") + d)
            return []
        else:
            d = sp(spaces)
            if is_R1:
                d += sp(5)
            d += sp(part.seq_len, "?")
            self._add_line(self._make_prefix("") + d)
            return []

    def _make_adapter_line(self, part, adapter, label):
        d = label
        d += sp(self.prefix_len + part.match_index + part.match_len + 1 - len(d))
        if (part == self.pair.r2):
            d += sp(5)
        d += (adapter[:part._rtrim - (4 if part == self.pair.r2 else 0)] + "..")
        self._add_line(d)

    def _make_part_errors(self, part):
        d = sp(part.match_index)
        errors = sp(part.seq_len)
        for e in part.match_errors:
            errors = errors[:e] + "!" + errors[e+1:]
        d += errors
        d += sp(self.target.n - len(d))
        if part.adapter_errors:
            d += " "
            if (part == self.pair.r2):
                d += sp(5)
                errors = sp(part._rtrim - 4)
            else:
                errors = sp(part._rtrim)
            for e in part.adapter_errors:
                errors = errors[:e] + "!" + errors[e+1:]
            if errors[0] == " ":
                errors = "|" + errors[1:]
            if errors[-1] == " ":
                errors = errors[:-1] + "|"
            d += errors
        self._add_line(self._make_prefix("-err!") + d)

    def _make_r1_rev(self):
        r1 = self.pair.r1
        d = "(revcomp)"
        d += sp(r1.match_index + self.prefix_len - len(d))
        d += r1.reverse_complement
        self._add_line(d)

    def _make_result(self, part):
        d = sp(part.match_index)
        l = "l={}".format(part.match_index)
        r = "r={}".format(part.match_index + part.match_len)
        leftover = part.match_len - len(l) - len(r) - 4
        if leftover < 0:
            l = l[2:]
            r = r[2:]
            leftover = max(0, part.match_len - len(l) - len(r) - 4)
        halfish = (leftover >> 1)
        bit = "^{}{}, {}{}^".format("-"*halfish,l,r,"-"*(leftover - halfish))
        d += bit
        self._add_line(sp(self.prefix_len) + d)

    def _add_line(self, line):
        if len(line) < self.max_len:
            line += sp(self.max_len - len(line))
        self.max_len = max(self.max_len, len(line))
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

    def _snip_lines(self, start, end):
        for i in range(len(self.lines)):
            line = self.lines[i]
            snip = line[start:end]
            replace = "  " if snip == sp(end-start) else ".."
            self.lines[i] = line[:start] + replace + line[end:]

    def make(self):
        self.lines = [ ]

        self._add_line("@" + self.pair.identifier)
        
        r2_bars = self._make_part(self.pair.r2)
        self.bars.append(r2_bars)
        if self.pair.r2.match_errors or self.pair.r2.adapter_errors:
            self._make_part_errors(self.pair.r2)
        if self.pair.r2.trimmed:
            self._make_adapter_line(self.pair.r2, reverse_complement(self.run.adapter_t), "RC(adapter_t)")

        r1_bars = self._make_part(self.pair.r1)
        self.bars.append(r1_bars)
        if self.pair.r1.match_errors or self.pair.r1.adapter_errors:
            self._make_part_errors(self.pair.r1)
        if self.pair.r1.trimmed:
            self._make_adapter_line(self.pair.r1, self.run.adapter_b, "adapter_b")

        if self.pair.mask and self.pair.r1.matched:
            self._add_line("")
            self._make_r1_rev()

        self._add_line("")
        self._add_line("")

        self._make_target_lines()

        if self.pair.r2.matched:
            self._make_result(self.pair.r2)
            self.bars.remove(r2_bars)
        if self.pair.r1.matched:
            self._make_result(self.pair.r1)
        else:
            self._add_line("")

        if self.pair.matched and self.pair.left > 100:
            self._snip_lines(self.prefix_len + 15, self.prefix_len + 85)

        self._make_summary()

        return "\n".join(self.lines)

def diagram(pair, run):
    return Diagram(pair, run).make()

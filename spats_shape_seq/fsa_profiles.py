from spats_shape_seq.profiles import TargetProfiles
from spats_shape_seq.spats import Spats
from spats_shape_seq.target import _Target

def strWithFileAtPath(path):
    with open(path, 'r') as f:
        return f.read()

def writeDataToFileAtPath(s, path):
    with open(path, 'wb') as f:
        f.write(s.encode())
    return True

# quick hack to compute approximate beta values from an .fsa file
# xref card #429
def compute(path):
    lines = strWithFileAtPath(path).split('\n')
    head = lines[0]
    lines = [l for l in lines[1:] if l]
    treated = []
    untreated = []
    for l in lines:
        bits = l.split('\t')
        areaRx, areaBg = float(bits[4]), float(bits[6])
        treated.append(int(areaRx))
        untreated.append(int(areaBg))
    treated_depths = []
    t = 0
    untreated_depths = []
    u = 0
    for i in range(len(treated)):
        t += treated[i]
        u += untreated[i]
        treated_depths.append(t)
        untreated_depths.append(u)
    spats = Spats()
    spats._run = spats.run
    tp = TargetProfiles(spats, _Target("x", "A"*(72 - 7), 0),
                        treated, untreated,
                        treated_depths, untreated_depths,
                        [], [])
    tp.compute()
    #print(tp.treated)
    #print(tp.untreated)
    #print(tp.treated_depths)
    #print(tp.untreated_depths)
    #print(tp.beta)
    out = [ head ]
    i = 0
    for l, b in zip(lines, tp.beta):
        out.append("{}\t{}".format(l, b))
    writeDataToFileAtPath('\n'.join(out), "{}.beta".format(path))


if __name__ == '__main__':
    import sys
    compute(sys.argv[1])

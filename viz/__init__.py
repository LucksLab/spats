
from viz.ui import SpatsViz

def run():
    s = SpatsViz()
    s.start()
    print("SpatsViz started...")
    s.waitFor()


import json
import os
import tempfile
import subprocess
import sys

import spats_shape_seq


class Plotter(object):

    def __init__(self):
        self._spats_path = os.path.normpath(os.path.join(os.path.dirname(spats_shape_seq.__file__), ".."))
        self.processes = []

    def stop(self):
        for proc in self.processes:
            if proc.returncode is None:
                proc.kill()

    def submit_plot(self, data, filename = ''):
        self.submit_plots([data], filename)

    def submit_plots(self, data, filename = ''):
        temp_file = tempfile.NamedTemporaryFile(delete = False)
        temp_file.write(json.dumps([data, filename]))
        proc = subprocess.Popen(["python", "viz/plotter.py", temp_file.name], cwd = self._spats_path)
        self.processes.append(proc)

    def _show_plot(self, figinfo):
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        plt.figure(1)

        idx = 0
        fig = figinfo[0]
        filename = figinfo[1]
        n = len(fig)
        for res in fig:

            idx += 1
            plt.subplot(int("{}1{}".format(n, idx)))

            for plot in res["data"]:
                plt.plot(plot["x"], plot["y"], plot["m"])

            if "xlim" in res:
                plt.xlim(res["xlim"])
            else:
                plt.xlim(0, max(plot["x"]))

            if "ylim" in res:
                plt.ylim(res["ylim"])

            plt.legend([ p.get("label", "") for p in res["data"] ])
            plt.title(res["type"])
            plt.xlabel(res["x_axis"])
            plt.ylabel(res["y_axis"])

        #def onclick(event):
        #    plt.close()
        #cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

        if filename:
            plt.gcf().canvas.get_default_filename = lambda: "spats_{}.png".format(filename)

        plt.show()

def show_plot(data_file):
    plot_data = json.loads(open(data_file, 'rb').read())
    Plotter()._show_plot(plot_data)

if __name__ == '__main__':
    show_plot(sys.argv[1])

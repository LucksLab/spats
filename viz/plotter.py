
import multiprocessing

class Plotter(object):

    def __init__(self):
        self.queue = multiprocessing.Queue()
        self.running = True
        self.fork()

    def fork(self):
        self.worker = multiprocessing.Process(target = self.handle_plots, args = [])
        self.worker.start()

    def handle_plots(self):
        try:
            while self.running:
                res = self.queue.get()
                if res:
                    worker = multiprocessing.Process(target = self._show_plot, args = (res,))
                    worker.start()
                else:
                    self.running = False
        except KeyboardInterrupt:
            pass

    def stop(self):
        self.queue.put(None)
        self.worker.join()

    def submit_plot(self, data, filename = ''):
        self.queue.put([[data], filename])

    def submit_plots(self, data, filename = ''):
        self.queue.put([data, filename])

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

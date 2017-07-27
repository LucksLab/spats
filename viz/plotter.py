
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

    def submit_plot(self, data):
        self.queue.put(data)

    def _show_plot(self, res):
        import matplotlib.pyplot as plt
        for plot in res["data"]:
            plt.plot(plot["x"], plot["y"], plot["m"])

        if "ylim" in res:
            plt.ylim(res["ylim"])

        plt.legend([ p.get("label", "") for p in res["data"] ])
        plt.title(res["type"])
        plt.xlabel(res["x_axis"])
        plt.ylabel(res["y_axis"])

        def onclick(event):
            plt.close()
        cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
        plt.show()

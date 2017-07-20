
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
        fig, ax = plt.subplots(figsize=(15, 10))
        plt.axes()
        plt.gca().add_patch(plt.Circle((42, 80), radius=0.02, fc = 'g', alpha = 0.2))
        plt.axis('scaled')
        plt.xlim([0,1])
        plt.ylim([0,1])

        def onclick(event):
            plt.close(fig)

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()

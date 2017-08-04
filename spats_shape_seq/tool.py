
import ConfigParser
import datetime
import multiprocessing
import os
import subprocess
import sys
import time

import spats_shape_seq
from spats_shape_seq import Spats
from spats_shape_seq.reads import ReadsData, ReadsAnalyzer
from viz.ui import SpatsViz


class SpatsTool(object):

    def __init__(self, path):
        self.path = path or os.getcwd()
        self.config = None
        self.cotrans = False
        self._skip_log = False
        self._no_config_required_commands = [ "viz" ]
        self._temp_files = []
        self._r1 = None
        self._r2 = None
        self._parse_config()

    def _parse_config(self):
        config_path = os.path.join(self.path, "spats.config")
        if not os.path.exists(config_path):
            return
        parser = ConfigParser.SafeConfigParser()
        parser.read(config_path)
        config = {}
        for section in parser.sections():
            if section not in config:
                config[section] = {}
            for name, value in parser.items(section):
                config[section][name] = value
        if config and 'spats' in config:
            self.config = config['spats']
            if 'cotrans' in self.config:
                cotrans = self.config['cotrans']
                self.cotrans = False if (cotrans == 'False') else bool(cotrans)

    def _spats_path(self):
        return os.path.normpath(os.path.join(os.path.dirname(spats_shape_seq.__file__), ".."))

    def _native_tool(self, tool_name):
        bin_path = os.path.join(self._spats_path(), "native", "bin", tool_name)
        return bin_path if os.path.exists(bin_path) else None

    def _add_note(self, note):
        self._notes.append(note)
        print(":{}".format(note))

    def _load_r1_r2(self, ret_r1):
        def decomp(rx):
            base, ext = os.path.splitext(rx)
            if ext.lower() == '.gz':
                self._add_note("decompressing {}".format(rx))
                out = base + ".tmp"
                subprocess.check_call("gzip -d -c {} > {}".format(rx, out), cwd = self.path, shell = True)
                self._temp_files.append(out)
                return out
            return rx
        self._r1, self._r2 = decomp(self.config['r1']), decomp(self.config['r2'])
        return self._r1 if ret_r1 else self._r2

    @property
    def r1(self):
        return self._r1 or self._load_r1_r2(True)

    @property
    def r2(self):
        return self._r2 or self._load_r1_r2(False)


    def _run(self, command):

        if not self.config and command not in self._no_config_required_commands:
            print("Missing spats.config")
            return

        self._notes = []

        start = time.time()
        hdlr = getattr(self, command, None)
        if not hdlr:
            print("Invalid command: {}".format(command))
            return
        hdlr()
        delta = time.time() - start

        if not self._skip_log:
            self._log(command, delta)

        for f in self._temp_files:
            if os.path.exists(f):
                os.remove(f)

    def _log(self, command, delta):
        stamp = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')
        with open(os.path.join(self.path, 'spats.log'), 'ab') as outfile:
            outfile.write("{} : {}, {}s\n".format(stamp, command, int(delta)))
            for note in self._notes:
                outfile.write("   - {}\n".format(note))
            outfile.write("\n")
                
    def reads(self):

        db_basename = 'reads.spats'
        db_name = os.path.join(self.path, db_basename)
        if os.path.exists(db_name):
            self._add_note("** removing previous {}".format(db_basename))
            os.remove(db_name)

        native_tool = self._native_tool('reads')
        if native_tool:
            self._add_note("using native reads")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, db_name], cwd = self.path)

        data = ReadsData('reads.spats')
        if not native_tool:
            self._add_note("using python reads")
            data.parse(self.config['target'], self.r1, self.r2)

        analyzer = ReadsAnalyzer(data, cotrans = self.cotrans)
        analyzer.process_tags()
        self._add_note("tags processed to {}".format(db_basename))

    def process(self):

        run_basename = 'run.spats'
        run_name = os.path.join(self.path, run_basename)
        if os.path.exists(run_name):
            self._add_note("** removing previous {}".format(run_basename))
            os.remove(run_name)

        native_tool = self._native_tool('cotrans')
        if self.cotrans and native_tool:
            self._add_note("using native cotrans processor")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, run_name], cwd = self.path)

        else:
            self._add_note("using python cotrans processor")
            spats = Spats(cotrans = self.cotrans)
            spats.addTargets(self.config['target'])
            spats.process_pair_data(self.r1, self.r2)
            spats.store(run_name)
        self._add_note("wrote output to {}".format(run_basename))

    def validate(self):
        run_basename = 'run.spats'
        run_name = os.path.join(self.path, run_basename)
        if not os.path.exists(run_name):
            print("Missing: {}".format(run_basename))
            return
        spats = Spats()
        spats.load(run_name)
        if spats.validate_results(self.r1, self.r2):
            self._add_note("Validation pass")
        else:
            self._add_note("Validation FAILURE")

    def viz(self):
        self._skip_log = True
        if sys.platform != "darwin":
            print("Invalid platform for viz UI: {}".format(sys.platform))
            return
        subprocess.check_call(["make", "vizprep"], cwd = self._spats_path())
        def viz_handler():
            sv = SpatsViz()
            sv.stop_on_disconnect = True
            sv.start()
            sv.waitFor()
        def uiclient_handler():
            bin_path = os.path.join(self._spats_path(), "bin", "UIClient.app", "Contents", "MacOS", "UIClient")
            try:
                subprocess.call([bin_path], cwd = self._spats_path())
            except KeyboardInterrupt:
                pass
        viz_worker = multiprocessing.Process(target = viz_handler, args = [])
        viz_worker.start()
        uiclient_worker = multiprocessing.Process(target = uiclient_handler, args = [])
        time.sleep(0.1)
        uiclient_worker.start()
        try:
            while viz_worker.is_alive() and uiclient_worker.is_alive():
                viz_worker.join(0.1)
                uiclient_worker.join(0.1)
        except KeyboardInterrupt:
            pass
        if viz_worker.is_alive():
            viz_worker.join()
        if uiclient_worker.is_alive():
            uiclient_worker.terminate()

        
        

def run(command, path = None):
    SpatsTool(path)._run(command)


if __name__ == '__main__':
    run(sys.argv[1])

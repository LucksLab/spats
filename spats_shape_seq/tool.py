
import ast
import ConfigParser
import csv
import datetime
import multiprocessing
import os
import subprocess
import sys
import time

import spats_shape_seq
from spats_shape_seq import Spats
from spats_shape_seq.reads import ReadsData, ReadsAnalyzer


class SpatsTool(object):

    def __init__(self, path):
        self.path = path or os.getcwd()
        self.config = None
        self.cotrans = False
        self._skip_log = False
        self._no_config_required_commands = [ "viz", "help", "doc" ]
        self._private_commands = [ "doc", "viz" ]
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
                self.cotrans = bool(ast.literal_eval(cotrans))

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
            base, ext = os.path.splitext(os.path.basename(rx))
            if ext.lower() == '.gz':
                self._add_note("decompressing {}".format(rx))
                out = os.path.join(self.path, base + ".tmp")
                subprocess.check_call('gzip -d -c "{}" > "{}"'.format(rx, out), cwd = self.path, shell = True)
                self._temp_files.append(out)
                self._sentinel("decompress {}".format('R1' if ret_r1 else 'R2'))
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


    def _run(self, args):
        if not args:
            print("Command required. Try 'spats_tool help'.")
            return

        command = args[0]
        self._command_args = args[1:]
        if not self.config and command not in self._no_config_required_commands:
            print("Missing spats.config")
            return

        self._notes = []

        self.start = time.time()
        hdlr = getattr(self, command, None)
        if not hdlr or (command in self._private_commands and spats_shape_seq._PUBLIC_RELEASE):
            print("Invalid command: {}".format(command))
            return
        try:
            hdlr()
            failure = False
        except Exception, e:
            print("** Command {} failed. ({})".format(command, e))
            failure = True

        if not failure and not self._skip_log:
            delta = self._sentinel("{} complete".format(command))
            self._log(command, delta)

        for f in self._temp_files:
            if os.path.exists(f):
                os.remove(f)

    def _sentinel(self, label):
        delta = time.time() - self.start
        self._add_note("{} @ {:.2f}s".format(label, delta))
        return delta

    def _log(self, command, delta):
        stamp = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')
        with open(os.path.join(self.path, 'spats.log'), 'ab') as outfile:
            outfile.write("{} : {}, {:.2f}s\n".format(stamp, command, delta))
            for note in self._notes:
                outfile.write("   - {}\n".format(note))
            outfile.write("\n")
                
    def reads(self):
        """Perform reads analysis on the r1/r2 fragment pairs, for use with the visualization tool.
        """

        db_name = self._reads_file()
        if os.path.exists(db_name):
            self._add_note("** removing previous reads file")
            os.remove(db_name)

        native_tool = self._native_tool('reads')
        if native_tool:
            self._add_note("using native reads")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, db_name], cwd = self.path)

        data = ReadsData(db_name)
        if not native_tool:
            self._add_note("using python reads")
            data.parse(self.config['target'], self.r1, self.r2)

        analyzer = ReadsAnalyzer(data, cotrans = self.cotrans)
        self._update_run_config(analyzer.run)
        analyzer.process_tags()
        self._add_note("tags processed to {}".format(os.path.basename(db_name)))

    def run(self):
        """Process the SPATS data for the configured target(s) and r1/r2 fragment pairs.
        """
        run_name = self._run_file()
        if os.path.exists(run_name):
            self._add_note("** removing previous run file")
            os.remove(run_name)

        native_tool = self._native_tool('cotrans')
        if self.cotrans and native_tool:
            self._add_note("using native cotrans processor")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, run_name], cwd = self.path)

        else:
            if native_tool:
                self._add_note("skipping native tool due to non-cotrans run")
            self._add_note("using python processor")
            spats = Spats(cotrans = self.cotrans)
            self._update_run_config(spats.run)
            spats.addTargets(self.config['target'])
            spats.process_pair_data(self.r1, self.r2)
            spats.store(run_name)
        self._add_note("wrote output to {}".format(os.path.basename(run_name)))

    def _update_run_config(self, run):
        for key, value in self.config.iteritems():
            if key in [ "r1", "r2", "target", "cotrans" ]:
                continue
            if hasattr(run, key):
                try:
                    val = ast.literal_eval(value)
                except:
                    val = value
                setattr(run, key, val)
                self._add_note("config set {} = {}".format(key, val))
            else:
                self._add_note("warning: unknown config {}".format(key))

    def _spats_file(self, base_name):
        return os.path.join(self.path, '{}.spats'.format(base_name))

    def _run_file(self):
        return self._spats_file('run')

    def _reads_file(self):
        return self._spats_file('reads')

    def validate(self):
        """Validate the results of a previous 'process' run against a second (slower) algorithm.
        """

        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be performed before validating")

        spats = Spats()
        spats.load(run_name)
        if spats.validate_results(self.r1, self.r2):
            self._add_note("Validation pass")
        else:
            self._add_note("Validation FAILURE")

    def viz(self):
        """Launch the visualization tool UI.
        """

        self._skip_log = True
        if sys.platform != "darwin":
            raise Exception("Invalid platform for viz UI: {}".format(sys.platform))

        subprocess.check_call(["make", "vizprep"], cwd = self._spats_path())
        def viz_handler():
            from viz.ui import SpatsViz
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

    def dump(self):
        """Dump data. Provide 'reads' or 'run' as an argument to dump the
        indicated type of data.
        """
        self._skip_log = True
        if not self._command_args:
            raise Exception("Dump requires a type (either 'reads' or 'run').")
        dump_type = self._command_args[0]
        handler = getattr(self, "_dump_" + dump_type, None)
        if not handler:
            raise Exception("Invalid dump type: {}".format(dump_type))
        handler()

    def _dump_reads(self):
        reads_name = self._reads_file()
        if not os.path.exists(reads_name):
            raise Exception("Reads must be run before attempting dump")
        data = ReadsData(reads_name)
        db = data.pair_db
        counts = db.tag_counts(1)
        total = float(db.count()) / 100.0
        keys = sorted(counts.keys(), key = lambda x : counts[x], reverse = True)
        data = [ [ key, float(counts[key])/total, counts[key] ] for key in keys ]
        output_path = os.path.join(self.path, 'reads.csv')
        self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], data)

    def _dump_run(self):
        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be run before attempting dump")

        spats = Spats()
        spats.load(run_name)
        profiles = spats.compute_profiles()
        headers = [ "L", "site", "nt", "f+", "f-", "beta", "theta", "rho", "c" ]
        data = []

        if self.cotrans:
            tgt = spats.targets.targets[0]
            tseq = tgt.seq
            for key in profiles.cotrans_keys():
                end = int(key.split('_')[-1])
                prof = profiles.profilesForTargetAndEnd(tgt.name, end)
                for i in range(end + 1):
                    data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.beta[i], prof.theta[i], prof.rho[i], prof.c ])
            output_path = os.path.join(self.path, '{}.csv'.format(tgt.name))
            self._write_csv(output_path, headers, data)
        else:
            for tgt in spats.targets.targets:
                tseq = tgt.seq
                end = len(tgt.seq)
                prof = profiles.profilesForTarget(tgt)
                data = []
                for i in range(end + 1):
                    data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.beta[i], prof.theta[i], prof.rho[i], prof.c ])
                output_path = os.path.join(self.path, '{}.csv'.format(tgt.name))
                self._write_csv(output_path, headers, data)


    def _write_csv(self, output_path, headers, data):
        with open(output_path, 'wb') as out_file:
            writer = csv.writer(out_file)
            if headers:
                writer.writerow(headers)
            for row in data:
                writer.writerow(row)
        self._add_note("Data dumped to {}".format(os.path.basename(output_path)))

    def doc(self):
        """Show the spats documentation.
        """

        self._skip_log = True
        subprocess.check_call(["make", "showdocs"], cwd = self._spats_path())

    def help(self):
        """Display available commands.
        """

        self._skip_log = True
        print("\nspats_tool v{}\nCommands:\n".format(spats_shape_seq._VERSION))
        for key in sorted(SpatsTool.__dict__.keys()):
            if key.startswith('_'):
                continue
            if key in self._private_commands and spats_shape_seq._PUBLIC_RELEASE:
                continue
            value = getattr(SpatsTool, key)
            if value.__doc__:
                print("  {}:  {}".format(key, value.__doc__.rstrip()))
        print("\nData commands require 'spats.config' to be present in the current working directory,")
        print("which requires a [spats] section header and should at least specify the 'target', 'r1',")
        print("and 'r2' configuration keys.\n")



def run(args, path = None):
    SpatsTool(path)._run(args)


if __name__ == '__main__':
    run(sys.argv[1:])

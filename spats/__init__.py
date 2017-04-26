
from config import spats_config
from spats import Spats
from util import reverse_complement

# helper function for simple cases
def spats(target_path, r1_path, r2_path, output_path, masks = ['RRRY', 'YYYR' ]):
     spats = Spats()
     spats.addTargets(target_path)
     spats.addMasks(*masks)
     spats.process_pair_data(r1_path, r2_path)
     spats.compute_profiles()
     spats.write_reactivities(output_path)


from functools import cache, cached_property
import json
from logging import getLogger
import os
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd

from . import path
from .sect import Section
from .seq import seq_to_unicode_array, write_fasta

logger = getLogger(__name__)

# Set environment variable DATAPATH to RNAstructure data tables.
data_path_key = "DATAPATH"
if os.environ.get(data_path_key) is None:
    rs_dir = "rnastructure"
    dt_dir = "data_tables"
    env_dir = os.path.dirname(os.path.dirname(sys.executable))
    for root, dirs, files in os.walk(env_dir):  # FIXME: this might not be portable
        if (os.path.basename(root) == dt_dir
                and os.path.basename(os.path.dirname(root)) == rs_dir):
            # RNAstructure data tables were found in root.
            # Set DATAPATH environment variable to root.
            os.environ[data_path_key] = root
            break
    else:
        # RNAstructure data tables were not found.
        # raise FileNotFoundError(os.path.join(rs_dir, dt_dir))
        pass


def run_command(cmd):
    import subprocess
    print(cmd)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8')


class RNAstructure(object):
    def __init__(self,
                 rnastructure_use_temp=None,
                 rnastucture_fold_args=None,
                 rnastructure_use_dms=None,
                 rnastructure_dms_min_unpaired_value=None,
                 rnastructure_dms_max_paired_value=None,
                 rnastructure_deltag_ensemble=None,
                 rnastructure_probability=None,
                 verbose=False,
                 ) -> None:
        self.rnastructure_path = rnastructure_path
        self.rnastructure_use_temp = rnastructure_use_temp
        self.rnastucture_fold_args = rnastucture_fold_args
        self.rnastructure_use_dms = rnastructure_use_dms
        self.rnastructure_dms_min_unpaired_value = rnastructure_dms_min_unpaired_value
        self.rnastructure_dms_max_paired_value = rnastructure_dms_max_paired_value
        self.rnastructure_deltag_ensemble = rnastructure_deltag_ensemble
        self.rnastructure_probability = rnastructure_probability
        self.verbose = verbose

        if self.rnastructure_use_dms is not None:
            assert self.rnastructure_dms_min_unpaired_value is not None, "If you want to use DMS, you need to specify the min and max values for the DMS data"
            assert self.rnastructure_dms_max_paired_value is not None, "If you want to use DMS, you need to specify the min and max values for the DMS data"

        if len(self.rnastructure_path) > 0:
            if self.rnastructure_path[-1] != '/':
                self.rnastructure_path += '/'
        self.directory = temp
        if os.path.exists(os.path.join(self.directory, 'ledger.json')):
            self.ledger = json.load(open(os.path.join(self.directory, 'ledger.json'), 'r'))
        else:
            self.ledger = {}
        os.makedirs(self.directory, exist_ok=True)

    def __make_files(self, temp_prefix='temp'):
        self.pfs_file = os.path.join(self.directory, temp_prefix + '.pfs')
        self.ct_file = os.path.join(self.directory, temp_prefix + '.ct')
        self.dot_file = os.path.join(self.directory, temp_prefix + '.dot')
        self.fasta_file = os.path.join(self.directory, temp_prefix + '.fasta')
        self.prob_file = os.path.join(self.directory, temp_prefix + '_prob.txt')

    def fit(self, sequence, reference='reference'):
        self.sequence = sequence
        self.__make_temp_folder()
        self.__make_files()
        self.__create_fasta_file(reference, sequence)

    def predict_ensemble_energy(self):
        cmd = f"{self.rnastructure_path}EnsembleEnergy {self.fasta_file} --sequence"
        splitted_output = run_command(cmd).split(' ')
        return float(splitted_output[splitted_output.index(f"kcal/mol\n\nEnsemble") - 1])

    def predict_partition(self, temperature_k=None):
        cmd = f"{self.rnastructure_path}partition {self.fasta_file} {self.pfs_file}"
        if temperature_k != None:
            cmd += ' --temperature ' + str(temperature_k)
        run_command(cmd)
        run_command(self.rnastructure_path + 'ProbabilityPlot ' + self.pfs_file + ' -t ' + self.prob_file)
        with open(self.prob_file, "r") as f:
            lines = f.readlines()
            out = {'i': [], 'j': [], 'p': []}
            for x in range(len(lines)):
                if x > 1:
                    ls = lines[x].split("\t")
                    out["i"] += [int(ls[0])]
                    out["j"] += [int(ls[1])]
                    out["p"] += [float(ls[2])]
        return self.__cast_pairing_prob(out)

    def predict_reference_deltaG(self):
        cmd = f"{self.rnastructure_path}Fold {self.fasta_file} {self.ct_file}"
        run_command(cmd)
        assert os.path.exists(
            self.ct_file), f"{self.ct_file} does not exist, check that RNAstructure works. Command: {cmd}"
        assert os.path.getsize(
            self.ct_file) != 0, f"{self.ct_file} is empty, check that RNAstructure works. Command: {cmd}"
        return self.__extract_deltaG_struct()

    def __make_temp_folder(self):
        isExist = os.path.exists(self.directory)
        if not isExist:
            os.makedirs(self.directory)
        return self.directory

    def __create_fasta_file(self, reference, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>' + reference + '\n' + sequence)
        temp_fasta.close()

    # cast the temp file into a dot_bracket structure and extract the attributes
    def __extract_deltaG_struct(self):
        run_command(f"{self.rnastructure_path}ct2dot {self.ct_file} 1 {self.dot_file}")
        temp_dot = open(self.dot_file, 'r')
        first_line = temp_dot.readline().split()
        # If only dots in the structure, no deltaG 

        if len(first_line) == 4:
            _, _, deltaG, _ = first_line
            deltaG = float(deltaG)
        if len(first_line) == 1:
            deltaG, _ = 0.0, first_line[0][1:]

        sequence = temp_dot.readline()[:-1]  # Remove the \n
        structure = temp_dot.readline()[:-1]  # Remove the \n
        return deltaG, structure

    def dump_ledger(self):
        json.dump(self.ledger, open(os.path.join(self.directory, 'ledger.json'), 'w'))

    def run(self, sequence):
        if sequence in self.ledger:
            deltaG, structure = self.ledger[sequence]['deltaG'], self.ledger[sequence]['structure']
        else:
            self.fit(sequence)
            deltaG, structure = self.predict_reference_deltaG()
            self.ledger[sequence] = {'deltaG': deltaG, 'structure': structure}
            for file in os.listdir(self.directory):
                if file.startswith('temp'):
                    os.remove(os.path.join(self.directory, file))
        if deltaG == 'void':
            deltaG = 0.0
        else:
            deltaG = float(deltaG)
        return {'deltaG': deltaG, 'structure': structure}


class RnaSection(object):
    """ RNA sequence or section thereof. """

    def __init__(self, title: str, section: Section):
        self.title = title
        self.section = section

    @property
    def ref(self):
        return self.section.ref

    @property
    def sect(self):
        return self.section.name

    @cached_property
    def seq(self):
        """ Sequence in RNA format. """
        return self.section.seq.tr()

    @property
    def seq_record(self):
        return self.section.ref_sect, self.seq

    @cached_property
    def index_from_one(self):
        """ Return a numerical index of the section starting from 1. """
        return pd.Index(self.section.positions - (self.section.end5 - 1))


class RnaProfile(RnaSection):
    """ Mutational profile of an RNA from a specific sample. """

    def __init__(self, title: str, section: Section, sample: str,
                 reacts: pd.Series):
        super().__init__(title, section)
        self.sample = sample
        self.reacts = reacts.reindex(self.section.positions)

    def get_ceiling(self, quantile: float) -> float:
        """ Compute the maximum reactivity given a quantile. """
        if not 0. < quantile <= 1.:
            raise ValueError("Quantile for reactivity ceiling must be in "
                             f"(0, 1], but got {quantile}")
        return np.nanquantile(self.reacts.values, quantile)

    def normalize(self, quantile: float):
        """ Normalize the reactivities to a given quantile. """
        ceiling = self.get_ceiling(quantile)
        if not 0. < ceiling <= 1.:
            raise ValueError("Reactivity ceiling must be in (0, 1], "
                             f"but got {ceiling}")
        return self.reacts / ceiling

    def winsorize(self, quantile: float):
        """ Normalize and winsorize the reactivities. """
        return pd.Series(np.clip(self.normalize(quantile), 0., 1.),
                         index=self.reacts.index)

    @cache
    def get_dir(self, out_dir: Path):
        """ Get the directory in which this RNA's files will go. """
        return path.builddir(path.ModSeg, path.SampSeg,
                             path.RefSeg, path.SectSeg,
                             top=out_dir, module=path.MOD_FOLD,
                             sample=self.sample, ref=self.ref, sect=self.sect)

    @cache
    def get_file(self, out_dir: Path, segment: path.Segment, **kwargs):
        """ Get the path to a file of the RNA sequence. """
        return self.get_dir(out_dir).joinpath(segment.build(**kwargs))

    def get_fasta(self, out_dir: Path):
        """ Get the path to the FASTA file of the RNA sequence. """
        return self.get_file(out_dir, path.FastaSeg,
                             ref=self.title, ext=path.FASTA_EXTS[0])

    def get_ct(self, out_dir: Path):
        """ Get the path to the FASTA file of the RNA sequence. """
        return self.get_file(out_dir, path.ConnectTableSeg,
                             struct=self.title, ext=path.CT_EXT)

    def get_dot(self, out_dir: Path):
        """ Get the path to the FASTA file of the RNA sequence. """
        return self.get_file(out_dir, path.DotBracketSeg,
                             struct=self.title, ext=path.DOT_EXTS[0])

    def get_dms(self, out_dir: Path):
        """ Get the path to the DMS data file of the RNA sequence. """
        return self.get_file(out_dir, path.DmsReactsSeg,
                             reacts=self.title, ext=path.DMS_EXT)

    def to_fasta(self, out_dir: Path):
        """ Write the RNA sequence to a FASTA file. """
        fasta = self.get_fasta(out_dir)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_dms(self, out_dir: Path, quantile: float):
        """ Write the DMS reactivities to a DMS file. """
        # The DMS reactivities must be numbered starting from 1 at the
        # beginning of the section, even if the section does not start
        # at 1. Renumber the section from 1.
        dms = pd.Series(self.winsorize(quantile), self.index_from_one)
        # Drop bases with missing data to make RNAstructure ignore them.
        dms.dropna(inplace=True)
        # Write the DMS reactivities to the DMS file.
        dms_file = self.get_dms(out_dir)
        dms.to_csv(dms_file, sep="\t", header=False)
        return dms_file


class Rna2dStructure(RnaSection):
    """ RNA secondary structure. """

    IDX_FIELD = "n"
    BASE_FIELD = "Base"
    PRIOR_FIELD = "n-1"
    NEXT_FIELD = "n+1"
    PAIR_FIELD = "Pair"
    POS_FIELD = "Position"

    def __init__(self,
                 title: str,
                 section: Section,
                 pairs: Iterable[tuple[int, int]]):
        super().__init__(section)
        self.title = title
        self.pairs = set(pairs)

    @property
    def header(self):
        return f"{self.title} {self.section.length}"

    @cached_property
    def partners(self):
        """ Return a Series of every position in the section and the
        position to which it pairs, or 0 if it does not pair. """
        # Initialize the series of pairs to 0 for every position.
        partners = pd.Series(0, index=pd.Index(self.section.positions,
                                               name=self.POS_FIELD))

        def add_pair(at: int, to: int):
            """ Add a base pair at position `at` to position `to`. """
            # Find the current pairing partner at this position.
            try:
                to2 = partners.loc[at]
            except KeyError:
                raise ValueError(f"Position {at} is not in {self.section}")
            if to2 == 0:
                # There is no pairing partner at this position: add it.
                partners.loc[at] = to
            elif to2 == to:
                # A previous partner matches the current partner.
                logger.warning(f"Pair {at}-{to} was given twice")
            else:
                # A previous partner conflicts with the current partner.
                raise ValueError(f"Position {at} was given pairs with both "
                                 f"{to2} and {to}")

        # Add all base pairs (in both directions) to the table.
        for pos1, pos2 in self.pairs:
            add_pair(pos1, pos2)
            add_pair(pos2, pos1)
        return partners

    @property
    def ct_data(self):
        """ Return the connectivity table as a DataFrame. """
        # Make an index the same length as the section and starting
        # from 1 (CT files must start at index 1).
        index = pd.RangeIndex(1, self.section.length + 1, name=self.IDX_FIELD)
        # Adjust the numbers of the paired bases (i.e. pairs > 0) such
        # that they also index from 1.
        pairs = self.partners.values.copy()
        pairs[pairs > 0] -= self.section.end5 - 1
        # Generate the data for the connectivity table.
        data = {
            self.BASE_FIELD: seq_to_unicode_array(self.seq),
            self.PRIOR_FIELD: index.values - 1,
            self.NEXT_FIELD: index.values + 1,
            self.PAIR_FIELD: pairs,
            self.POS_FIELD: self.section.positions,
        }
        # Assemble the data into a DataFrame.
        return pd.DataFrame.from_dict({field: pd.Series(values, index=index)
                                       for field, values in data.items()})

    @property
    def ct_text(self):
        """ Return the connectivity table as text. """
        data = self.ct_data.reset_index()
        return f"{self.header}\n{data.to_string(index=False, header=False)}\n"


class RnaState(Rna2dStructure):
    """ RNA sequence and mutation rates. """

    def __init__(self,
                 title: str,
                 section: Section,
                 pairs: Iterable[tuple[int, int]],
                 sample: str,
                 mus: pd.Series):
        """
        Parameters
        ----------
        mus: pd.DataFrame
            Mutation rates of the RNA molecule.
        """
        super().__init__(title, section, pairs)
        self.sample = sample
        self.mus = mus


class RnaEnsemble(RnaSection):
    pass



def fold(ct_file):
    pass

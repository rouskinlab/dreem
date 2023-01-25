from dreem.draw import manipulator, util, plotter
import pandas as pd
import numpy as np
from dreem.util.dump import sort_dict, flatten_json

class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'])
    """

    attr_list = ['name','samples']

    def __init__(self, data=None, samples=None, min_cov_bases=0, filter_by='sample') -> None:
        """Creates a Study object.

        Args:
            data (dict or list[dict], optional): A dictionary or list of dictionaries containing the data to be loaded. One dictionary corresponds to a fastq file. Defaults to None.
            samples (List[str], optional): List of samples to load. Defaults to None.
            min_cov_bases (int, optional): Minimum number of base coverage for a row to be filtered-in. Defaults to 0.
            filter_by (str, optional): Filter rows by sample or study. When filtered by study, if a row passes the filter, rows with the same 'construct', 'section' and 'cluster' fields for all other samples have a sufficient base coverage. Defaults to 'sample'.            

        Example:
            >>> study = Study(data = {'sample':'mysample',{'construct1': {'section1': {'cluster1': {'mut_bases': [100], 'cov_bases': [1000]}}}}},
                              samples=['A1', 'B2', 'B3'], 
                              min_cov_bases=1000, 
                              filter_by='study')
        """
        self.samples = samples
        if data is not None:
            self.df = pd.DataFrame()
            if type(data) is not list:
                data = [data]
            print('Turning data into a dataframe...')
            for sample in data:
                print(sample, end='... ')
                self.df = pd.concat([self.df, pd.DataFrame(flatten_json(sort_dict(sample)))], axis=0)
            print('Done.')
            print('Setting dataframe...')
            self.set_df(self.df, min_cov_bases=min_cov_bases, filter_by=filter_by, samples=samples)
            print('Done.')
        else:
            self.df = None

    @classmethod
    def from_dict(cls, di:dict):
        """Set attributes of this Study object from a dictionary.

        Args:
            di (dict): a dictionary containing keys such as ['name','description','samples','label','conditions'].

        Returns:
            Study: a study object.

        Example:
        >>> di = {'name':'temperature','samples':['A1','B2','B3']}
        >>> study = Study().from_dict(di)
        >>> print(study.name, study.samples)
        temperature ['A1', 'B2', 'B3']
        """
        for attr in cls.attr_list:
            try: 
                di[attr]
            except: 
                di[attr]=None
        return cls(di['name'], di['samples'])

    def set_df(self, df, min_cov_bases=0, filter_by='sample', samples=None):
        df.reset_index(inplace=True, drop=True)
        self.df = df
        
        if not 'worst_cov_bases' in self.df.columns:
            self.df['worst_cov_bases'] = self.df['cov_bases'].apply(lambda x: min(x))
        
        self.df = self.df[self.df['worst_cov_bases'] >= min_cov_bases]
        if samples is not None:
            assert type(samples) is list, 'samples must be a list'
            self.df = self.df[self.df['sample'].isin(samples)]
        
        if filter_by == 'study':
            self.filter_by_study(inplace=True)
        
        for attr in ['sample','construct']:
            self.df[attr] = self.df[attr].astype(str)
        
        for attr in ['section','cluster']:
            if attr not in self.df.columns:
                self.df[attr] = 0
        
        self.df['deltaG'] = self.df['deltaG'].apply(lambda x: 0.0 if x == 'void' else float(x))
    
    
    def filter_by_study(self, inplace=False):
        df = self.df.groupby(['construct', 'section', 'cluster']).filter(lambda x: len(self.df['sample'].unique()) == len(x['sample'].unique()))
        if inplace:
            self.df = df
        return df.copy()

    def get_df(self, **kwargs):
        return manipulator.get_df(self.df, **kwargs)

    def get_samples(self):
        return self.df['sample'].unique()

    def get_constructs(self, sample:str):
        return self.df[self.df['sample'] == sample]['construct'].unique()

    def get_sections(self, sample:str, construct:str):
        return self.df[(self.df['sample'] == sample) & (self.df['construct'] == construct)]['section'].unique()

    def get_clusters(self, sample:str, construct:str, section:str):
        return self.df[(self.df['sample'] == sample) & (self.df['construct'] == construct)& (self.df['section'] == section)]['cluster'].unique()
       
    def load_studies(studies_file_path:str):
        return load_studies(studies_file_path)



    def mutation_fraction(self, **kwargs)->dict:
        """Plot the mutation rates as histograms.

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
            base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
            base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by expected base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            show_ci(bool, optional): Show confidence intervals. Defaults to True.
            savefile(str, optional): Path to save the plot. Defaults to None.
            use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
            title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.

        Returns:
            dict: Figure and data of the output plot.
        """

        return plotter.mutation_fraction(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k in list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)}), **{k:v for k,v in kwargs.items() if k in plotter.mutation_fraction.__code__.co_varnames})

    def deltaG_vs_mut_rates(self, **kwargs)->dict:
        """Plot the mutation rate of each paired-expected base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
            min_cov_bases (int, optional): Filter rows by a minimum threshold for base coverage. Defaults to 0.
            base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
            base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by expected base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            savefile(str, optional): Path to save the plot. Defaults to None.
            use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
            title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            **kwargs: Additional arguments to pass to filter rows by. Ex: flank='flank_1' will keep only rows with flank=flank_1. 

        Returns:
            dict: Figure and data of the output plot.
        """
        return plotter.deltaG_vs_mut_rates(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k in list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)}), **{k:v for k,v in kwargs.items() if k in plotter.deltaG_vs_mut_rates.__code__.co_varnames})

    
    def exp_variable_across_samples(self, **kwargs)->dict:
        """Plot the mutation rate of each paired-expected base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            construct (str): Construct of your row.
            experimental_variable (str): x axis column value, must be a per-sample attribute.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (str): Cluster of your row.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            savefile (str, optional): Path to save the plot. Defaults to None.
            use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
            title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.

        Returns:
            dict: Figure and data of the output plot.
        """
        return plotter.exp_variable_across_samples(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k in  list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)}), **{k:v for k,v in kwargs.items() if k in plotter.exp_variable_across_samples.__code__.co_varnames})
    
    def auc(self, **kwargs)->dict:
        """Plot the AUC for each mutation profile of the selected data. 

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
            min_cov_bases (int, optional): Filter rows by a minimum threshold for base coverage. Defaults to 0.
            base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
            base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by expected base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            savefile(str, optional): Path to save the plot. Defaults to None.
            use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
            title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.
            **kwargs: Additional arguments to pass to filter rows by. Ex: flank='flank_1' will keep only rows with flank=flank_1. 
    
        Returns:
            dict: Figure and data of the output plot.
        """
        return plotter.auc(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k in list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)}), **{k:v for k,v in kwargs.items() if k in plotter.auc.__code__.co_varnames})


    def mutation_fraction_delta(self, **kwargs)->dict:
        """Plot the mutation rate difference between two mutation profiles.
        sample1: sample of the first mutation profile.
        sample2: sample of the second mutation profile.
        construct1: construct of the first mutation profile.
        construct2: construct of the second mutation profile.
        section1: section of the first mutation profile.
        section2: section of the second mutation profile.
        cluster1: cluster of the first mutation profile.
        cluster2: cluster of the second mutation profile.
        base_index1: base index of the first mutation profile.
        base_index2: base index of the second mutation profile.
        base_type1: base type of the first mutation profile.
        base_type2: base type of the second mutation profile.
        base_pairing1: base pairing of the first mutation profile.
        base_pairing2: base pairing of the second mutation profile.
        savefile(str, optional): Path to save the plot. Defaults to None.
        use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
        title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.        
        """
        
        df1 = manipulator.get_df(self.df, **{k[:-1]:v for k,v in kwargs.items() if k.endswith('1') and k[:-1] in list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)})
        assert len(df1)>0, 'No rows found for the first mutation profile.'
        assert len(df1)==1, 'More than one row found for the first mutation profile.'
        df2 = manipulator.get_df(self.df, **{k[:-1]:v for k,v in kwargs.items() if k.endswith('2') and k[:-1] in list(self.df.columns)+ list(manipulator.get_df.__code__.co_varnames)})
        assert len(df2)>0, 'No rows found for the second mutation profile.'
        assert len(df2)==1, 'Only one row should be selected for the second mutation profile.'
        return plotter.mutation_fraction_delta(pd.concat([df1, df2]).reset_index(drop=True), **{k:v for k,v in kwargs.items() if k in plotter.mutation_fraction_delta.__code__.co_varnames})

    def base_coverage(self, **kwargs):
        """Plot the base coverage of several constructs in a sample.

        Args:
            samp (str): Sample of your rows.
            constructs (List[str]): Constructs of your rows.
            section (str): Region of your row.
            cluster (int, optional): Cluster of your row. Defaults to 0. 
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``. Can be a series of 0-indexes (ex: [43,44,45,48]), 'roi', 'all', or a unique sequence (ex: 'ATTAC')
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            base_paired (bool, optional): Base-pairing predicition to plot. Defaults to None.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            show_ci (bool, optional): Show confidence interval on the histogram. Defaults to True.
            savefile (str, optional): Path to save the plot. Defaults to None.
            use_iplot(bool, optional): Use iplot instead of plot (for Jupyter notebooks). Defaults to True.
            title(str, optional): Title of the plot. Defaults to None, in which case a standard name is given.

        Returns:
            dict: Figure and data of the output plot.

        """
        return 0# plotter.base_coverage(self._df, **kwargs)



def load_studies(studies_file_path:str)->dict[str:Study]:
    """Read formatted file with samples, and turn it into a dataframe containing studies.

    Args:
        studies_file_path (str): path+title of the csv file containing the samples.

    Returns:
        (pd.DataFrame): studies of the csv file, indexed by study.
    """

    studies_dict, studies_data = {}, pd.read_csv(studies_file_path)

    for col in studies_data.groupby('name')[Study.attr_list]:
        solo_item = lambda x: x[0] if len(set(x)) == 1 else x  
        studies_dict[col[0]] = {attr: solo_item(list(col[1][attr])) for attr in (Study.attr_list)} 

    return {k:Study.from_dict(v) for k,v in studies_dict.items()}

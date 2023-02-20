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

    def __init__(self, data=None, min_cov_bases=0, filter_by='sample') -> None:
        """Creates a Study object.

        Args:
            data (dict or list[dict] or pandas.DataFrame, optional): Data to use. Can be a dictionary or list of dictionaries containing DREEM-output jsons, or directly a pandas dataframe. Defaults to None.
            min_cov_bases (int, optional): Minimum number of base coverage for a row to be filtered-in. Defaults to 0.
            filter_by (str, optional): Filter rows by sample or study. When filtered by study, if a row passes the filter, rows with the same 'reference', 'section' and 'cluster' fields for all other samples have a sufficient base coverage. Defaults to 'sample'.            

        Example:
            >>> study = Study(data = {'sample':'mysample',{'reference1': {'section1': {'cluster1': {'mut_bases': [100], 'cov_bases': [1000]}}}}},
                              min_cov_bases=1000, 
                              filter_by='sample')
        """
        if data is not None:
            
            df = pd.DataFrame()              
            
            # If data is a list of json, concatenate them into a single dataframe
            if type(data) is not pd.DataFrame:
                print('Turning data into a dataframe...')
                
                if type(data) is not list:
                    data = [data]
                    
                for sample in data:
                    print(sample['sample'], end='... ')
                    df = pd.concat([df, pd.DataFrame(flatten_json(sort_dict(sample)))], axis=0)
                
                print('Done.')
            
            # Use the dataframe (loaded or created from json)
            print('Setting dataframe...')
            self.set_df(df, min_cov_bases=min_cov_bases, filter_by=filter_by)
            print('Done.')
            
        else:
            self.df = None

    
    def set_df(self, df, min_cov_bases=0, filter_by='sample'):
        
        self.df = df.reset_index(drop=True)
                
        self.df = self.df[self.df['worst_cov_bases'] >= min_cov_bases]
        
        if filter_by == 'study':
            self.filter_by_study()
        
        for attr in ['sample','reference']:
            self.df[attr] = self.df[attr].astype(str)
        
        for attr in ['section','cluster']:
            if attr not in self.df.columns:
                self.df[attr] = 0
        
        self.df['deltaG'] = self.df['deltaG'].apply(lambda x: 0.0 if x == 'void' else float(x))
    
    
    def filter_by_study(self):
        df = self.df.groupby(['reference', 'section', 'cluster']).filter(lambda x: len(self.df['sample'].unique()) == len(x['sample'].unique()))
        self.df = df

    def get_df(self, **kwargs):
        return manipulator.get_df(self.df, **kwargs)

    def get_samples(self):
        return self.df['sample'].unique()

    def get_references(self, sample:str):
        return self.df[self.df['sample'] == sample]['reference'].unique()

    def get_sections(self, sample:str, reference:str):
        return self.df[(self.df['sample'] == sample) & (self.df['reference'] == reference)]['section'].unique()

    def get_clusters(self, sample:str, reference:str, section:str):
        return self.df[(self.df['sample'] == sample) & (self.df['reference'] == reference)& (self.df['section'] == section)]['cluster'].unique()
       

    def mutation_fraction(self, **kwargs)->dict:
        """Plot the mutation rates as histograms.

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            reference (list, int, str, optional): Filter rows by reference (list of references or just a reference). Defaults to None.
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
        """Plot the mutation rate of each paired-expected base of the ROI for each reference of a sample, w.r.t the deltaG estimation.

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            reference (list, int, str, optional): Filter rows by reference (list of references or just a reference). Defaults to None.
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
        """Plot the mutation rate of each paired-expected base of the ROI for each reference of a sample, w.r.t the deltaG estimation.

        Args:
            reference (str): reference of your row.
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
            reference (list, int, str, optional): Filter rows by reference (list of references or just a reference). Defaults to None.
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
        reference1: reference of the first mutation profile.
        reference2: reference of the second mutation profile.
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
        """Plot the base coverage of several references in a sample.

        Args:
            sample (str): Sample of your rows.
            reference (List[str]): references of your rows.
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


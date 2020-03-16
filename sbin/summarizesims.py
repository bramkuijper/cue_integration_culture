#!/usr/bin/env python3

import os
import re
import sys
import pandas as pd
import multiprocessing as mp

# class to summarize simulation output
# where csv files with data + parameters are parsed
# so that a summary file is generated with the last
# line of the data + all parameters + filename for each file
class SummarizeSims:
    def __init__(self
                 ,path="."
                 ,pattern=r"(sim.*\d+$|iter.*\d+)"
                 ,recursive=False
                 ,testing=False
                 ,sep=";"
                 ,n_process=3
                 ,parameters_first=False
                 ,posthoc_function=None
                 ,max_number_files=None
                 ):
        """
        Initializes the summarize sims class

        Parameters
        ----------
        path : str, optional
            The path in which files need to found. The default is ".".
        pattern : str, optional
            Regular expression matching the filenames that need to be 
            processed. The default is r"sim|iter".
        recursive : boolean, optional
            Whether the directory should be searched recursively. The
            default is False
        sep : str, optional
            The separator used in the file. The default is ";"
        n_process : int
            The number of nodes working on this (multiprocessing)
        parameters_first : boolean, optional
            Whether the listing of the parameters occurs 
            at the start of the file
        posthoc_function: function, optional
            Function to be called after processing data of each simulation

        Returns
        -------
        None.

        """
        
        self.first = True
        self.init = False
        self.path = path
        self.num_semicol_header = 0
        self.pattern = pattern
        self.recursive = recursive
        self.sep = sep
        self.testing = testing
        self.n_process = n_process
        
        self.posthoc_function = posthoc_function
        self.full_data = None


        if self.testing:
            print("traversing path '" + self.path + "'")

        file_ctr = 0

        # list of files where work is needed
        self.file_list = []
        
        # loop through directory tree to get a list of files 
        # which match the regex
        for root, dirname, files in os.walk(self.path):
            
            for file in files:
                if re.search(self.pattern,file) is not None:
                    self.file_list.append(os.path.join(root, file))

        # truncate if only a number of files needed
        self.file_list = self.file_list[:max_number_files]

        # now perform multiprocessing to get a dataframe
        # with all records
        # always leave one cpu process empty
        if self.n_process >= os.cpu_count() - 1:
            self.n_process = os.cpu_count() - 1

        # make pool object for multiprocessing
        pool = mp.Pool(processes=self.n_process) 

        result = pool.map(
                self.analyze_file
                ,self.file_list
                )

        self.full_data = pd.concat(result)

        if self.full_data is not None:
            
            self.output_postprocess()
            self.output_full_data()

    def output_postprocess(self):
        """
        Post collection data cleaning: e.g., strip whitespace
        from column names

        Returns
        -------
        None.

        """
        self.full_data.rename(str.strip
                       ,inplace=True
                       ,axis="columns")
        
        

    def output_full_data(self):
        
        

        # write the data to stdout
        if self.full_data.shape[0] > 0:
            print(self.full_data.to_csv(
                    path_or_buf=None
                    ,sep=";"
                    ,index=False))

    # analyze parameters at the end of the file
    def analyze_parameters_old(self, lines):
        """
        Runs through the parameter part and analyzes it

        Parameters
        ----------
        lines : list
            list of strings, where each string is a line from the original
            simulation file.
        
        Returns
        -------
        None.

        """
    
        pars = {}
    
        for line in lines:
            mobj = line.split(self.sep)
    
            if len(mobj) > 1:
                pars[mobj[0]] = mobj[1]
    
        return(pars)

    def analyze_data(
            self
            ,filename
            ,linenumber_start=0
            ,linenumber_end=None):

        the_data = pd.read_csv(
                filepath_or_buffer=filename
                ,skiprows=linenumber_start
                ,nrows=linenumber_end
                ,sep=self.sep)
    
        return(the_data.iloc[[-1]])
    
    def analyze_parameters_new(
            self 
            ,linenumber_start
            ,filename
            ,linenumber_end=None):

        assert(type(linenumber_start) == type(30))
        assert(type(filename) == type("adf"))
        
        if linenumber_end != None:
            assert(linenumber_end > linenumber_start)

        # read in the parameter data as if it were
        # a csv file
        param_data = pd.read_csv(
                filepath_or_buffer=filename
                ,sep=self.sep
                ,skiprows=linenumber_start
                ,nrows=linenumber_end
                ,header=None
                ,dtype=str
                ,usecols=[0,1]).T

        param_data.columns = param_data.iloc[0]
        param_data.drop(index=0,axis=0, inplace=True)

        return(param_data) 
    
    # processes the first line headers
    # when making line headers for initial values
    def process_first_line(self, line):
    
        # get the column names and split them into a list
        line_cols = line.strip().split(self.sep)
    
        new_cols = ""
    
        for colname in line_cols:
            if not colname:
                continue
    
            new_cols += colname + "_t_0" + self.sep
    
        return(new_cols)

    def find_last_data_line(self, filename):

        with open(filename) as infile:
            linelist = infile.readlines()

            # make a reverse range and loop through lines
            linerange = range(len(linelist)-1,0,-1)

            for line_idx in linerange:
                if re.search(r"^\d",linelist[line_idx]) is not None:
                    return(line_idx)
            else:
                raise ValueError("Cannot find any data in file " + filename)

    
    def analyze_file(self,filename):

        # two options (ignoring empty lines)
        # 1. header line, data, then parameters
        # 2. parameters, header line, data 

        # parameters recognizable as they contain
        # non-numerical data at start. Hence multiple
        # non-numerical starts consecutively means parameters

        non_data_lines_start = []
        data_start = None

        # open the file
        with open(filename) as infile:

            # loop through lines from start
            for line_idx, line in enumerate(infile):
                # skip empty lines
                if line in ["\n","\r\n"]:
                    continue

                # check whether line contains parameters
                # or column headings
                if re.search(r"^\s*\d",line) is None:
                    non_data_lines_start.append(line_idx)
                else:
                    # ok we hit a line containing numerical
                    # data
                    # we can break as we know enough
                    data_start = line_idx
                    break
            else:
                # ok no data was found, as a break should
                # have occurred somewhere
                # we will ignore this file
                return None


        # we are done sniffing the data, we can now make a choice
        # if more than one non_data_line, we got a parameter listing
        # at the start
        if len(non_data_lines_start) > 1:

            # parameter listing at the start, call analyze parameters function
            parameters = self.analyze_parameters_new(
                    linenumber_start=non_data_lines_start[0]
                    ,linenumber_end=non_data_lines_start[-2] # not the last element as that is the data header
                    ,filename = filename)
            
            data_posthoc = None

            if self.posthoc_function != None:
                data_posthoc = self.posthoc_function(filename)

            data = self.analyze_data(
                    filename=filename
                    ,linenumber_start=non_data_lines_start[-1]
                    )
          

        else: # parameters at the end apparently

            # get a list of all non-data-lines at the end of the file
            last_data_line = self.find_last_data_line(filename=filename)
            
            # collect the data
            data = self.analyze_data(
                    linenumber_start=non_data_lines_start[-1]
                    ,linenumber_end=last_data_line
                    ,filename=filename)

            data_posthoc = None

            if self.posthoc_function != None:
                data_posthoc = self.posthoc_function(filename)
            
            parameters = self.analyze_parameters_new(
                    linenumber_start=last_data_line + 1
                    ,linenumber_end=None
                    ,filename = filename)
            
                   
        # combine both parameters and data
        data_params_combined = pd.concat(
                [parameters.reset_index(drop=True)
                    ,data.reset_index(drop=True)
                    ,data_posthoc.reset_index(drop=True)
                    ]
                ,axis=1)
            
        # add current filename
        data_params_combined["file"] = filename

        return(data_params_combined)
        

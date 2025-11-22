import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.gridspec import GridSpec
from natsort import natsorted
from datetime import datetime
import os
from .TAhelpers import metadata_cond, regions_dict, find_region, split_label_list, div_calc, dict_val_fraction, return_2std, find_qnt,raw2div


# Analyze the trimers and produce a scatter plot of trimers diversity per position in documents.
class TrimerAnalysis():
    def __init__(self, 
                 trimer_input: pd.DataFrame,
                 relabel: list = None,
                 ):
    
        """
        1. trimer_input : str -> input trimer csv file location or pd.DataFrame object
        2. relabel : list of tuples -> list of tuples to rearange the dataset labels in the format of [("index_name1", index_int1), ("index_name2", index_int2), ...]
                                    for example: [("subject_id",2), ("ab_target",0)]
        3. trimer_count : str -> source of trimer count in positions and reporiores. 'diversity' will only include highly frequent trimers and 'richness' will take
                                into consideration all of the occurnces.
        """
    
        # Importing the trimers output
        if isinstance(trimer_input, (pd.DataFrame)):
            self.raw_data = trimer_input
        
        else:
            try:
                self.raw_data = pd.read_csv(trimer_input, index_col=0)
            except:
                raise SyntaxError("Invalid trimer input path")
            
        # Chaning the J region pos_aa format to "j.X" from "Xj" for better presentation
        self.raw_data.loc[self.raw_data.region == "j", "first_aa"] = self.raw_data[self.raw_data.region == "j"].first_aa.apply(lambda X : X[-1] + "." + X[:-1])
        
        # Can re-label the dataset based on defined list of toples
        if relabel is not None:
            labels_splited = self.raw_data.label.str.split(".")
            index_rearanged = labels_splited.map(lambda X : X[relabel[0][1]])

            if len(relabel) > 1:
                for index_i in relabel[1:]:
                    index_rearanged = index_rearanged + "." + labels_splited.map(lambda X : X[index_i[1]])
            
            self.raw_data.label = index_rearanged
        
    # Creating dataframe of occurnces of trimers, divrasity and richness in each position
    def diversity_analysis(self,
                          save_csv : bool = False,
                          div_orders: list = [0, 1]) -> pd.DataFrame:

        """
        save_csv : bool -> save the output csv data
        div_orders : list -> Orders of diversity to calculate.
        """
        self.div_orders = div_orders

        # Creating diversity calculation list
        div_list = []
        div_columns = [f"div{i}" for i in div_orders]
        div_unq_columns = [f"div{i}_unq" for i in div_orders]
        div_unq_grouped = [[f"div{i}", f"div{i}_unq"] for i in div_orders]
        for i in div_orders:
            div_list += [f"div{i}", f"div{i}_unq"]

        # Grouping the data by label>pos>trimer and summing the number of counts of each trimer
        raw_grouped = self.raw_data.groupby(["label", "first_aa","trimer"])["ncount"].sum().reset_index()

        column_names = ["label", "aa_frac"] + div_list # name of the columns
        template_df = pd.DataFrame(index=natsorted(self.raw_data.first_aa.unique()), columns=column_names).reset_index(names="pos")
        divresity_dflst = []

        # Itirating over the different documents
        for i_label in raw_grouped.label.unique():
            temp_grouped = raw_grouped[raw_grouped.label == i_label].iloc[:,1:]
            label_template = template_df.copy()

            # Itiriating over the different positions
            for i_pos in label_template.pos:
                temp_pos = temp_grouped[temp_grouped.first_aa == i_pos]

                # Incase of 0 diversity == 1 (log(0)=1) so we put 0 and nan in advance
                if len(temp_pos) == 0:
                    label_template.loc[label_template.pos == i_pos, "label"] = i_label
                    label_template.loc[label_template.pos == i_pos, div_columns] = 0
                    label_template.loc[label_template.pos == i_pos,["aa_frac"] + div_unq_columns] = np.nan

                else:
                    fraction_series = temp_pos[["trimer", "ncount"]].set_index("trimer", drop=True)["ncount"]
                    fraction_series = (fraction_series/fraction_series.sum()).sort_values(ascending=False)
                    label_template.loc[label_template.pos == i_pos, ["label", "aa_frac"]] = i_label, {i:round(float(j),3) for i,j in zip(fraction_series.index, fraction_series.values)}

                    for i,j in zip(div_orders, div_unq_grouped):
                        div_i = div_calc(array=fraction_series.values.tolist(), order=i)
                        label_template.loc[label_template.pos == i_pos, j] = div_i, np.array(fraction_series[:div_i].index)

            # appeding the docuemnt dataframe to the list of dataframes
            divresity_dflst.append(label_template)
        
        self.output_diversity = pd.concat(divresity_dflst)
        self.output_diversity["region"] = self.output_diversity["pos"].apply(find_region) 

        if save_csv:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_results\\trimers_diversity_{}.csv'.format(time)

            if os.path.exists("trimers_results") is False:
                    os.mkdir("trimers_results")

            self.output_diversity.to_csv(save_path)
        
        return self.output_diversity
    

    # Counting the number of repotiores and position each trimer exists in
    # Counting the number of repotiores and position each trimer exists in
    def trimer_presence(self,
                        index_labels : dict,
                        drop_cdr3j : bool = True, 
                        div_order : list = None,
                        diversity : list = [1],
                        save_csv : bool = False
                        ) -> pd.DataFrame:
            
            """
            index_labels : dict -> Dictionary with the index and metadata names that construct the label.
            drop_cdr3j : bool -> to drop cdr3 and j region from the analysis output.
            div_order : str -> Source of trimer calculation count (diversity order, need to be in div_orders argument of TrimerAnalysis).
            diversity : list-like -> array of int > 0 values, order of diversity to calculate (presence of amino acid in said diveristy).
            save_csv : bool -> Save the outout csv.
            """

            try:
                itrimer_source = {"richness":"div0_unq", "diversity":"div_unq"}
                input_source = f"div{div_order}_unq"

            except:
                raise Exception(f"invalid 'div_order' input, {div_order} not in {self.div_orders}.") 
            
            trimers_info = self.raw_data
            if drop_cdr3j:
                trimers_info = trimers_info[trimers_info.region.isin(['cdr1', 'fw2', 'cdr2', 'fw3'])]
            
            template_index = trimers_info.trimer.unique()
            # template_cnames = ["label", "trimer", "n_pos", "unique_pos"]
            template_dataframe = pd.DataFrame(index=template_index)
            
            final_list = []
            for i_label in trimers_info.label.unique():
                # grouping by -> reseting index -> fixing index levels -> setting 'trimers' columns as index -> deleting the index name
                temp_trimers = trimers_info[trimers_info.label == i_label].groupby(["label", "trimer"]).agg({"first_aa":["nunique", "unique", dict_val_fraction], "ncount":"sum"}).reset_index().droplevel(level=1, axis=1).set_index("trimer").rename_axis(None)
                temp_trimers.columns = ["label", "div0_npos", "div0_unqpos", "div0_pcount", "div0_ncount"]
                temp_trimers = template_dataframe.join(other = temp_trimers, how="left")
                temp_trimers.label = i_label
                            
                # Appending the i_label grouped datagrame to list
                final_list.append(temp_trimers)
                 
            self.trimer_count_final = pd.concat(final_list, axis=0)
            self.trimer_count_final.loc[self.trimer_count_final.div0_npos.isnull(), "div0_npos"] = 0
            self.trimer_count_final.div0_npos = self.trimer_count_final.div0_npos.astype("int")
            self.trimer_count_final.reset_index(drop=False, names="trimer", inplace=True)

            if isinstance(diversity,  list):
                #try:
               for i in diversity:
                        temp_div = trimers_info.groupby(["label","first_aa"]).apply(raw2div, div=i, include_groups=False).reset_index()
                        test_gp = temp_div.groupby(["trimer","label"]).agg({"first_aa":["nunique", "unique", dict_val_fraction], "ncount":"sum"}).reset_index()
                        gp_cols = [f"div{i}_npos", f"div{i}_unqpos", f"div{i}_pcount", f"div{i}_ncount"] 
                        test_gp.columns = ["trimer", "label"] + gp_cols
                        test_gp.index = test_gp.trimer + "_" + test_gp.label

                        self.trimer_count_final.index = self.trimer_count_final.trimer + "_" + self.trimer_count_final.label
                        self.trimer_count_final = pd.merge(left=self.trimer_count_final, 
                                                           right=test_gp[gp_cols], 
                                                           left_index=True, 
                                                           right_index=True, 
                                                           how="left")
                        
            # Creating metadata columns according to our index_labels dict
            for ic in index_labels.keys():
                 self.trimer_count_final[index_labels[ic]] = self.trimer_count_final.label.str.split(".").apply(lambda X : X[ic])

            if save_csv:
                time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
                save_path = 'trimers_results\\trimers_presence_{}.csv'.format(time)

                if os.path.exists("trimers_results") is False:
                        os.mkdir("trimers_results")

                self.trimer_count_final.to_csv(save_path)

            return self.trimer_count_final.reset_index(drop=True)
    
    # plotting diversity/richness values across the amino acids positions
    def plot_positions(self,
                       label_index : int,
                       filter_datasets : int = 0,
                       #source : str = "diversity",
                       save_fig : bool = False,
                       sep : bool = True) -> plt.figure:
        
        """
        label_index : int -> Divide the data based on unique values of the label index (when the label is divided by ".").
        filter_datasets : int -> drop repotiores with less than 'filter_datasets' trimers 
        source : str -> source of the presentation 'diversity' or 'richness'
        save_fig : bool -> save the plot.
        sep : bool -> seperate 'j' and 'cdr3' region from the rest in the presentation.
        """
        
        #try:
        #    itrimer_source = {"richness":"rich", "diversity":"div"}
        #    trimer_source = itrimer_source[source]
        #except:
        #    raise Exception("invalid 'source' argument input, need to be string 'richness' or 'diversity'") 
        
        for div_order in self.div_orders:
            trimer_source = f"div{div_order}"
            print(f"Plotting Position Analysis of Diversity of Order {div_order}:")

            # Creating a filter for low values datasets
            labels_sums = self.raw_data.groupby("label").agg({"ncount":"sum"}).reset_index().sort_values(by="ncount")
            labels_ok = labels_sums[labels_sums["ncount"] > filter_datasets].label.values
            self.labels_small = labels_sums[labels_sums["ncount"] <= filter_datasets].label.values
            print(f"> Dropped small labels {filter_datasets} < : {self.labels_small}")

            dataset = self.output_diversity.copy()
            dataset = dataset[dataset["label"].isin(labels_ok)]

            # Creating colormapping dict for different labels 
            dataset["handle"] = dataset.loc[:,"label"].str.split(".").map(lambda X : X[label_index])
            colors = [cm.to_hex(plt.cm.tab10(i)) for i in range(len(dataset.handle.unique()))]
            chandles_dic = {i:j for i,j in zip(dataset.handle.unique(), colors)}

            #fig, axs = plt.subplots(2, 2, figsize=(25,10))
            fig = plt.figure(figsize=(25,10))
            fig_shape = (2,2)

            if sep:
                dic_datasets = {0:dataset[(dataset.region != "cdr3") & (dataset.region != "j")],
                                1:dataset[dataset.region == "cdr3"],
                                2:dataset[dataset.region == "j"]}
                
                ax1 = plt.subplot2grid(shape=fig_shape, loc=(0,0), rowspan=1, colspan=2)
                ax2 = plt.subplot2grid(shape=fig_shape, loc=(1,0), rowspan=1, colspan=1)
                ax3 = plt.subplot2grid(shape=fig_shape, loc=(1,1), rowspan=1, colspan=1)
                ax_list = [ax1, ax2, ax3]

                ax1_title = "CDR3 & J Excluded"

            else:
                ax1 = plt.subplot2grid(shape=fig_shape, loc=(0,0), rowspan=2, colspan=2)
                ax_list = [ax1]
                dic_datasets = {0:dataset}
                ax1_title = "All Regions"
            
            for idt, axi, ittl in zip(dic_datasets, ax_list, [ax1_title, "CDR3 Region", "J Region"]):
                dataset = dic_datasets[idt]
                ax_xticks = dataset["pos"].values

                for j in dataset.label.unique():
                    temp_handle = j.split(".")[label_index]
                    temp_df = dataset[dataset.label == j]

                    axi.scatter(temp_df["pos"], 
                                temp_df[trimer_source],
                                label=temp_handle,
                                color=chandles_dic[temp_handle],
                                s=10,
                                alpha=0.4)

                axi.tick_params(axis='x', labelrotation=90) # rotating xticks labels
                axi.set_xlabel("Amino Acid Position") 
                axi.set_ylabel(f"Trimer Diversity (Order {div_order})")
                axi.set_title(ittl)

                yaxis_ticks = axi.get_yticks()
                yaxis_interval = (abs(yaxis_ticks.max()) + abs(yaxis_ticks.min()))/(len(yaxis_ticks)-1)
                axi.set_yticks(ticks = np.arange(yaxis_ticks.min(), yaxis_ticks.max()+yaxis_interval*2, yaxis_interval)[1:-1])

                # Unique legends labels
                handles, labels = axi.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                axi.legend(by_label.values(), by_label.keys())

                # Changing xtick labels to red if the trimers is from CDR region
                def int_func(val):
                    try:
                        return int(val)
                    except:
                        return int(val.split(".")[1])
        
                if axi is ax1:
                    range_cdr = np.arange(27,39).tolist() + np.arange(56,66).tolist()
                    for xtick, xcolor in zip(ax1.get_xticklabels(), ["tab:red" if int_func(i) in range_cdr else "black" for i in ax_xticks]):
                        xtick.set_color(xcolor)
            
            if save_fig:
                time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
                save_path = 'trimers_figures\\trimers_div{}_positions_{}.png'.format(div_order, time)

                if os.path.exists("trimers_figures") is False:
                        os.mkdir("trimers_figures")

                fig.savefig(save_path, bbox_inches='tight')

            plt.tight_layout()
            plt.show()


# Presenting the number of trimers per position in different documents.
def vis_distribution(trimer:str,
                     drop_rows:int = 100,
                     region:list = None,
                     positions:list = None,
                     found_in_all:bool = True,
                     label_index:list = None,
                     pos_trimer:bool = False,
                     save_fig:bool = False):
    
    """
    trimer: str -> trimer source location, will be importaed as pd.DataFrame.
    drop_rows: int -> dropping sub-datasets with less than 'drop_rows' unique sequences.
    region: list -> Trimers that are found in the BCR heavy chain regions.
    positions: array-like -> Trimers that are found in the BCR heavy chain amino acid positions.
    found_in_all: bool -> Does the trimer found in all of the sub-datasets.
    label_index: tuple -> format of [(str, int),(,),...], plotting the sub-datasets according to the index label which divide the label
                          by '.' the [0] string will be used as legend title in the plot and the [1] int will decide on the string
                          division and the relvent sagment.
    pos_trimer: bool -> adding the position of the first aa to the trimer label.
    """

    # Importing the trimers dataset
    try:
        data = pd.read_csv(trimer, index_col=0)
        #data_uclones = pd.read_csv(unique_clones, index_col=0)

    except:
        raise Exception("Invalid data path.")
    
    # Adding amino acid position to the start of the trimer label
    if pos_trimer:
        data["trimer"] = data["first_aa"].astype("str") + "_" + data["trimer"] + "_" + data["region"]
    
    # Filtring by region and positions based on the 'region' and 'positions' function arguments
    low_vals = (data.groupby("label")["ncount"].max() > drop_rows) 
    data_filt = data.copy()[data.label.isin(low_vals[low_vals == True].index)]
    print(f"Dropped low count (<100) sub-datasets: {low_vals[low_vals==False].index}")

    cnd_dict = {}
    for cnd, cnd_str in zip([region, positions], ["region", "first_aa"]):
        if cnd is None:
            cnd_dict[cnd_str] = np.full((data_filt.shape[0],), fill_value=True, dtype=bool)
            
        elif cnd is not None:
            cnd_dict[cnd_str] = data[cnd_str].isin(cnd)
    
        data_filt = data_filt[cnd_dict[cnd_str]]

    # Creating a dataframe with all of the frequencies of the possible trimers for each sub-dataset (defined by 'label')
    data_filt = data_filt.groupby(["label", "trimer"])["ncount"].sum().reset_index()
    template_df = pd.DataFrame(index=data.trimer.unique())

    # Concantinaing the infromation about trimer frequency across the columns (horizently) - including missing values of all possible trimers
    first_i = True
    for i in data_filt.label.unique():
        temp_df = data_filt[data_filt.label == i].copy().reset_index(drop=True)
        temp_sum = temp_df.ncount.sum()
        temp_df[i] = temp_df.ncount / temp_sum
        temp_df = temp_df.sort_values(by=i, ascending=False).set_index("trimer")[i].to_frame(name=i)

        if first_i:
            result_df = pd.concat([template_df, temp_df], axis=1)
            first_i = False
        else:
            result_df = pd.concat([result_df, temp_df], axis=1)
    
    #setting the X axis (trimer) by avg fraction
    result_df = result_df.loc[result_df.median(axis=1).sort_values(ascending=False).index,:] 

    # dropping null rows (any-value) if argument found_in_all is True
    plot_input = result_df 
    if found_in_all:
        plot_input = result_df.dropna(axis=0, how="any")

    # number of trimers (size of x axis)
    n_trimers = len(plot_input.index) 
    
    # Incase we want to split the data labels by index
    column_names = plot_input.columns
    actual_label_index = plot_input.columns.str.split(".").map(len)[0] - 1
    if (label_index is not None):
        index_dict = {li[1]:column_names.str.split(".").map(lambda X : X[li[1]]).unique() for li in label_index}
        index_labels = ".".join([i[0] for i in label_index])

        uqlabels_df = metadata_cond(index_dict)
        unique_labels = [".".join(uqlabels_df.iloc[i,:].values) for i in range(0,uqlabels_df.shape[0])]
        
        #assiging colors to labels
        colors = [cm.to_hex(plt.cm.tab20(i)) for i in range(uqlabels_df.shape[0])]
        labels_dic = {i:j for i,j in zip(unique_labels, colors)}
        
        # renaming the index for the x-axis labels
        plot_input_old = plot_input.index
        if pos_trimer == True:
            plot_input.index= [i.split("_")[1] + f" ({i.split("_")[0]})" for i in plot_input.index]

        #for index_label in unique_labels:
        #    pass

    # Plotting the data
    fig, ax = plt.subplots(figsize = (40,7))
    for i_label in column_names:
        if label_index is None:
            temp_label = ""
            temp_color = "tab:blue"

        else:
            temp_label = None
            for i in label_index:
                if temp_label is None:
                    temp_label = f"{i_label.split(".")[i[1]]}"
                else:
                    temp_label += f".{i_label.split(".")[i[1]]}"
            
            temp_color = labels_dic[temp_label]

        ax.scatter(x = plot_input.index, 
                   y = plot_input[i_label],
                   color = temp_color,
                   label = temp_label,
                   s = 5,
                   alpha = 0.8)

    # Plot properties
    ax.axhline(y=0, color="black", alpha=0.5, lw=1, ls="--")
    ax.set_xlabel('Trimer')
    ax.tick_params(axis='x', labelrotation = 90, labelsize=7)
    ax.set_ylabel("Trimer Fraction")
    ax.grid(True, alpha=0.5, zorder=0)
    ax.set_xlim(-1, n_trimers)
    ax.text(x=0, y=plt.yticks()[0].max(), s=f"Number of trimers = {n_trimers}")

    # Changing xtick labels to red if the trimers is from CDR region
    if pos_trimer:
        range_cdr = ["cdr1", "cdr2", "cdr3"]
        for xtick, xcolor in zip(ax.get_xticklabels(), 
                                         ["tab:red" if i.split("_")[2] in range_cdr else "black" for i in plot_input_old]):
            xtick.set_color(xcolor)
        
    # Getting legend without duplicated values
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), title=index_labels)

    if save_fig:
        time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
        save_path = 'trimers_figures\\trimers_distribution_{}.png'.format(time)

        if os.path.exists("trimers_figures") is False:
                os.mkdir("trimers_figures")

        fig.savefig(save_path, bbox_inches='tight')

    return plot_input


# class that analyzes
class medians_info():
    def __init__(self,
                raw_trimers : pd.DataFrame,
                lpa_sig : pd.DataFrame,
                index_unique : dict,
                index_name : str):
                #trimers_diversity : pd.DataFrame,
                #trimers_array : np.array):
        """
        raw_trimers : pd.DataFrame -> raw trimers dataframe, created via TrimerCreate.Trimer.create().
        lpa_sig : pd.DataFrame -> lpa distances signautres of the trimers across the documents, created via the class object sig_df of the TrimerLPA() initiated class.
        index_unique : dict -> dictionary of the metadata index in the label and it's unique values. for example ({0:['Non-SpikeB', 'Spike+MemB']} for the first index whihch
                               contains the unique values of the antibody target unique values.
        index_name : str -> label of the metdata for example time-point, subject_id or maybe antibody type. should be in accordance with the dataset naming scheme
        """
        
        # Importing the index information and lpa signatures dataframe
        index = list(index_unique.keys())[0]
        self.unique_metadata = list(index_unique.values())[0] + ["joined"]
        self.index_name = index_name
        trimer_lpa = lpa_sig

        # creating template for the signatures median df
        self.trimer_sig_medians = pd.DataFrame(index = lpa_sig.index)

        # calcualting the trimer median per ab_target type (sn, sp)
        for i in self.unique_metadata:
            i_label = f"{i}_median"

            if i in self.unique_metadata[:-1]:
                i_cond = (trimer_lpa.columns.str.split(".").map(lambda X : X[index]) == i)
                i_series = trimer_lpa.loc[:,i_cond].median(axis=1)
                i_median = pd.DataFrame(index = i_series.index , data={i_label:i_series.values})

            else:
                i_median = pd.DataFrame(index =  trimer_lpa.median(axis=1).index , data={i_label: trimer_lpa.median(axis=1).values})

            i_median = i_median.sort_values(by=i_label, ascending=False)

            self.trimer_sig_medians = self.trimer_sig_medians.join(other=i_median, how="right")

        # Reseting index and chaning columns names
        self.trimer_sig_medians.reset_index(names="trimer", inplace=True)
        self.trimer_sig_medians.columns = ["trimer"] + ["_".join([i,j]) for i in self.unique_metadata for j in ["median"]]

        # Joining the signatures dataframe to the raw_trimers grouped
        # List of functions and columns for the groupby aggrgate function (dict vals -> in TAhelpers.py)
        agg_functions = ["nunique", "unique", dict_val_fraction] 
        in_columns = raw_trimers.columns.tolist()
        in_columns.remove("trimer")
        in_columns.remove("ncount")

        # dictionary of all the columns and appling the function of agg_function (line 6)
        agg_dict = {i:agg_functions for i in in_columns}
        agg_dict["ncount"]="sum"

        # Creating the dataset
        raw_grouped = raw_trimers.groupby("trimer").agg(agg_dict).reset_index()

        # dropping level of columns
        raw_grouped = raw_grouped.droplevel(level=1, axis=1)
        i_columns = ["trimer"] + ["_".join([i,j]) for i in raw_grouped.columns[1:-1].unique() for j in ["nunique","unique","frac"]]  + ["ncount_sum"]
        raw_grouped.columns = i_columns

        # joining the trimer mdian distances to the precentile dataframe
        self.raw_grouped_medians = raw_grouped.merge(right=self.trimer_sig_medians, left_on="trimer", right_on="trimer", how="left")

    # adding statistical information to the "on" median data, such as precentile or above/between or lower than 2std 
    def filter_sort(self,
                    on : str,
                    sort_ascending : bool = False,
                    filter_how : str = "qnt",
                    save_csv : bool = False) -> pd.DataFrame:
        """
        sort_by : str -> by what column to sort the dataframe.
        sort_ascending : str ->  how to sort the the dataframe, True for ascending and False for descending.
        filter_how : str -> how to filter, by 25% precentile (qnt) or 2 standard deviation (2std).
        save_csv : bool -> to save the csv output.
        """

        self.filt_column = on
        cond_sort_by = (isinstance(on, str) & (on in self.raw_grouped_medians.columns[-6:]))
        cond_sort_how = (isinstance(sort_ascending, bool))
        cond_filter_how = (isinstance(filter_how, (str, None)) & (filter_how in["qnt", "2std"]))
        
        cond_array = [cond_sort_by, cond_sort_how, cond_filter_how]
        if all(cond_array) is False:
            print(cond_array)
            raise TypeError("Invalid argument input, verify and re-run the function")
        

        self.medians_filtered = self.raw_grouped_medians.copy()
        if filter_how == "2std":
            value_2std =  2 * self.medians_filtered[on].std()
            self.medians_filtered["_".join([on, filter_how])] = self.medians_filtered[on].apply(return_2std, std=value_2std)
            
        
        if filter_how == "qnt":
            #quantiles = self.medians_filtered["Non-SpikeB_median"].describe()[3:].values
            quantiles =  self.medians_filtered[on].quantile([0, 0.25, 0.5, 0.75, 1])
            self.medians_filtered["_".join([on, filter_how])] = self.medians_filtered[on].apply(find_qnt, quantile=quantiles)

        self.medians_filtered.sort_values(by=on, ascending=sort_ascending, inplace=True)

        if save_csv:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_medians\\trimers_medians_{}.csv'.format(time)

            if os.path.exists("trimers_medians") is False:
                    os.mkdir("trimers_medians")
            
            self.medians_filtered.to_csv(save_path)
    
        return self.medians_filtered
    
    
    def plot_median_hist(self,
                         trimer_pres : pd.DataFrame,
                         div : int = 0,
                         xy_limit : tuple = None,
                         save_fig : bool = False):
    
        """
        trimer_pres : pd.DataFrame -> output dataframe of the TrimerAnalyze.trimer_presence() method.
        xy_limit : tuple -> limiting the x,y axis ticks to specific value in the format of (x_min, x_max, y_min, y_max)
        save_fig : bool -> bool True to save the figure or False to discard.
        """

        trimers_relv = trimer_pres[trimer_pres[f"div{div}_npos"].notnull()].trimer.unique()
        print(f"> {len(trimers_relv)} unique trimers in divresity of order {div}. Producing histogram.")
        
        # Getting the relecent metadata labels as defined in theclass initiation index_unique argument
        type_labels = self.unique_metadata
        len_tlabels = len(type_labels)
        
        # Defining figure object 
        fig, axs = plt.subplots(1, len_tlabels, figsize=(6*len_tlabels, 5))
        y_ticks = []
        x_ticks = []
        
        # plotting the figure with subplots
        for i,j in zip(axs, type_labels):
            i_vector = self.medians_filtered.loc[self.medians_filtered.trimer.isin(trimers_relv) ,f"{j}_median"].values
            sns.histplot(data=i_vector, bins=100, stat="percent", ax=i, kde=True)

            i.set_ylabel("Precentage %")
            i.set_xlabel(f"'{j}' Median Distance")
            i.tick_params(axis='x', labelrotation=45)
            i.set_title(f"'{j}' Domain Distances Medians Distributions (diversity {div})", size=10)
            
            for si in 2*[-i_vector.std(), i_vector.std()]:
                i.axvline(x=si, color="tab:red", alpha=0.5, lw=1, ls="--", zorder=10)
            

            x_ticks += [float(tick.get_text().replace("−","-")) for tick in i.get_xticklabels()]
            y_ticks += [float(tick.get_text().replace("−","-")) for tick in i.get_yticklabels()]
             
        # if xy_limit defined -> takes max and min for y and x axis in order to limit the subplots y and x ticks to the same limits
        try:
            x_ticks = [xy_limit[0], xy_limit[1]]
            y_ticks = [xy_limit[2], xy_limit[3]]

        except:
            x_ticks = [min(x_ticks), max(x_ticks)]
            y_ticks = [min(y_ticks), max(y_ticks)]
        
        for i in axs:
                i.set_xlim(x_ticks[0], x_ticks[1])
                i.set_ylim(y_ticks[0], y_ticks[1])

        #saving figure
        if save_fig:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_medians\\medians_hist_{}.png'.format(time)

            if os.path.exists("trimers_medians") is False:
                    os.mkdir("trimers_medians")

            fig.savefig(save_path, bbox_inches='tight')
        
        plt.show()

    def plot_median_presence(self,
                             trimer_pres : pd.DataFrame,
                             div_pres : int = 0,
                             stat : str = "qnt",
                             mask_zeros : bool = True,
                             save_fig : bool = False,
                             return_data : bool = True):
        """
        trimer_pres : pd.DataFrame -> output dataframe of the TrimerAnalyze.trimer_presence() method.
        div_pres : int -> what order to diversity to calculate, need to be analyzed in via the TrimerAnalyze.trimer_presence() method.
        mask_zeros : bool -> removing the non-present trimers from the data and plot.
        save_fig : bool -> bool True to save the figure or False to discard.
        return_data : bool -> bool to return 
        """

        if stat not in ["qnt", "std", "2std"]:
            print("> stat argument invalid input, performing the analysis for 2 time standatd diviation (std).")
            stat = "std"

        # labels of ab_targets
        type_labels = self.unique_metadata
        len_tlabels = len(type_labels)
        
        # setting universal y axis limit
        xy_medians_values = self.medians_filtered[[i + "_median" for i in self.unique_metadata + ["joined"]]]
        y_max = xy_medians_values.max().max()
        y_min = xy_medians_values.min().min()

        fig, axs = plt.subplots(1, len_tlabels, figsize=(6*len_tlabels, 5), constrained_layout=True)
        dict_top_trimers = {}

        for i,j in zip(axs, type_labels):
            # need:
            # 1. list of trimers per label
            # 2. medians of said trimers
            # 3. diviresity npos of said trimers
            # datasets: raw_trimers, trimers_presence and medians_filt

            label_name = self.index_name #class initiation argument
             # tdf_medians -> medians values of the sub-dataset
             # tdf_npos -> mean occurnces in positions of the trimers in wanted diveristy order and sub-dataset 
             # tdf_concat -> the dataset to plot [0] is median distance (y-axis) and [1] is the existance in mean positions in order of diversity

            trimer_count = trimer_pres

            zero_switch = {1: (trimer_count[f"div{div_pres}_npos"] != 0),
                           0: ([True] * trimer_count.shape[0])}

            if j in self.unique_metadata[:-1]:
                tdf_medians = self.medians_filtered[["trimer", f"{j}_median"]].set_index("trimer")
                tdf_npos = trimer_count.loc[(trimer_count[f"div{div_pres}_npos"] >= 1) & (trimer_count[label_name] == j),].groupby("trimer").agg({f"div{div_pres}_npos":"mean"}).dropna()
                tdf_concat = pd.concat([tdf_medians, tdf_npos], axis=1).loc[tdf_npos.index,:]

            else:
                 tdf_medians = self.medians_filtered[["trimer", f"{type_labels[-1]}_median"]].set_index("trimer")
                 tdf_npos = trimer_count.loc[(trimer_count[f"div{div_pres}_npos"] != 0)].groupby("trimer").agg({f"div{div_pres}_npos":"mean"}).dropna()
                 tdf_concat = pd.concat([tdf_medians, tdf_npos], axis=1)
            
            # information for the plot: precentiles, precentiles numbers, number of positions, colors for high medians and high trimer count
            tdf_describe = tdf_concat.iloc[:,0].describe()
            tdf_std = tdf_describe["std"]
            stats_dict = {"std":[pd.Series([-tdf_std, tdf_std]), "std"],
                          "2std":[pd.Series([2*-tdf_std, 2*tdf_std]), "2std"],
                          "qnt":[tdf_describe[4:-1], "Precentile"]}
        
            
            j_precentiles = stats_dict[stat][0]
            j_medians = tdf_concat.iloc[:,0].values
            j_npos = tdf_concat.iloc[:,1].values
            j_colors = []

            for trimer in j_medians:
                if trimer >= j_precentiles.values[-1]:
                    j_colors.append("tab:red")
                
                elif trimer <= j_precentiles.values[0]:
                    j_colors.append("tab:green")
                
                else:
                    j_colors.append("tab:blue")

            j_count = j_colors.count("tab:red")

            # saving the information into dictionary for further use (the method returns that dict)
            tdf_concat["color"] = j_colors
            output_df = tdf_concat.sort_values(f"{j}_median", ascending=False)
            dict_top_trimers[j] = output_df[output_df["color"] == "tab:red"]
                
            # plotting the dataset in the sub-plot
            i.scatter(x=j_npos ,y=j_medians, alpha=0.25, color=j_colors, label="absolute")
            i.set_xlabel("Mean Occurces in Amino Acid Positions")
            i.set_ylabel("Mean Trimer LPA Median Distance")
            i.set_title(f"'{j}' Amino Acid Median LPA Distance vs Amino Acid Occurnces", size=9)
            i.tick_params(axis='y', labelrotation=45)
            i.set_ylim(y_min, y_max)
            
            # Orgenizing y axis labels
            og_ticks = list(i.get_yticks())
            interval = abs(og_ticks[0]) - abs(og_ticks[1])
            new_ticks = [og_ticks[0] - interval] + og_ticks + [og_ticks[-1]+interval]
            i.set_ylim(new_ticks[1], new_ticks[-2])

            # Subplots text and vertical line for precentiles
            i.annotate(f"-- {stats_dict[stat][1]} line", color="tab:red", xy=(1, 0.9), xycoords='axes fraction', fontsize=10,
                xytext=(-5, 5), textcoords='offset points', ha='right', va='bottom')
            i.annotate(f"{j_count} top trimers", color="tab:red", xy=(1, 0.935), xycoords='axes fraction', fontsize=10,
                xytext=(-5, 5), textcoords='offset points', ha='right', va='bottom')
            
            # precentiles lines marking
            for pr,val in zip(j_precentiles.index, j_precentiles.values):
                i.axhline(y=val, zorder=2, color="tab:red", ls="--", alpha=0.5)

        #saving figure
        if save_fig:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_medians\\medians_presence_div{}_{}.png'.format(div_pres,time)

            if os.path.exists("trimers_medians") is False:
                os.mkdir("trimers_medians")

            fig.savefig(save_path, bbox_inches='tight')
        
        plt.show()
        print("> Non-precense trimers were dropped from the figure, the top trimers count != 25% of the trimers \n> Non-precense trimers still have median distance due to the scoring method.")
        
        if return_data:
            return dict_top_trimers
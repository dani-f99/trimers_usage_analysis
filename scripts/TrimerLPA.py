from sklearn.decomposition import PCA as skPCA
from sklearn.preprocessing import StandardScaler
from .LPA import Corpus, sockpuppet_distance


import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime


# Performs latent domain LPA analysis on the trimers frequencies whitin different immune repotiores, the input is the TrimerCreate.create_trimer() input.
class TrimerLPA():
    def __init__(self, 
                 input_data,
                 min_treshold:int = None,
                 keep_zero_freq:bool = True):
        
        """
        input_data: str -> variable of pd.DataFrame or csv file location of  trimmer output from the Trimers.create_trimers() method.
        min_treshold: int -> treshold for number of trimers in documents, lower values from this treshold will be dropped.
        keep_zero_freq: -> Keeping trimers with zero occurnces, the LPA algorithm drops them automaticly.
        """

        if isinstance(input_data, pd.DataFrame):
            raw_trimer = input_data
        else:
            raw_trimer = pd.read_csv(input_data, index_col=0)
        
        self.lpa_input = raw_trimer[["label","trimer","ncount"]]
        self.lpa_input.columns = ["document", "element", "frequency_in_document"]
        self.lpa_input = self.lpa_input.groupby(["document","element"]).sum().reset_index()

        if isinstance(min_treshold, int):
            doc_max = self.lpa_input.groupby("document").max().reset_index()
            self.low_docs = doc_max[doc_max.frequency_in_document < min_treshold].document.values
            self.high_docs = doc_max[doc_max.frequency_in_document >= min_treshold].document.values

            self.lpa_input = self.lpa_input[self.lpa_input.document.isin(self.high_docs)]
            print(f"> The following document havn't met the treshold ({str(min_treshold)}):")
            print(doc_max.loc[doc_max.document.isin(self.low_docs), ["document", "frequency_in_document"]])

        if keep_zero_freq:
            unique_trimers =  np.sort(self.lpa_input.element.unique())
            unique_labels = self.lpa_input.document.unique()
            template_df = pd.DataFrame(index=unique_trimers)

            temp_dfs = []
            for label in unique_labels:
                input_df = self.lpa_input[self.lpa_input.document == label].set_index("element")["frequency_in_document"]
                conc_df = pd.concat([template_df, input_df], axis=1).fillna(0).astype("int")
                conc_df["document"] = label
                temp_dfs.append(conc_df)

            self.lpa_input = pd.concat(temp_dfs, axis=0).reset_index(drop=False, names="element")[["document", "element", "frequency_in_document"]]

        # creating corpus and domain 
        corpus = Corpus(self.lpa_input)
        dvr = corpus.create_dvr()

        # sorting the domain by element value
        dvr_sort = dvr.sort_values(by="element").reset_index(drop=True)

        # defining epsilon for KLDe distance calculation
        epsilon_frac = 2
        epsilon = 1 / (len(dvr) * epsilon_frac)

        # creating distances signatures for each document
        signatures = corpus.create_signatures(epsilon=epsilon, sig_length=None, distance="KLDe")
        sig_list = [i.sort_index() for i in signatures]

        self.sig_df = dvr_sort.copy()["element"].to_frame()

        # Creating dataframe of signatures
        for i in sig_list:
           self.sig_df = self.sig_df.merge(right=i.to_frame(), how="left", left_on="element", right_index=True) 
            
        self.sig_df.set_index(keys="element", drop=True, inplace=True)
        self.sig_df = self.sig_df[np.sort(self.sig_df.columns)]


    # Heatmap plot: KLDe distances of each dataset from the shared domain
    def plot_heatmap(self,
                     trimers_origin:str,
                     mask_heatmap : bool = True,
                     save_fig : bool = True,
                     hide_heamapx : bool = True) -> plt.Figure:
        """
        trimers_origin : str -> origin of the trimers dataset, used for plot x axid label.
        mask_heatmap : bool -> to mask NaN values.
        save_fig : bool -> save the output figure.
        """
        
        domain_df = self.sig_df
        fhight = domain_df.shape[1]/5
        fig, ax = plt.subplots(1,1, figsize=(20,fhight))


        mask_df = np.full(domain_df.shape,False)
        if mask_heatmap:
            mask_df = domain_df.isnull()

        
        sns.heatmap(data=domain_df.T,
        cbar_kws={'label': 'KLDe Distance'},
        cmap="turbo",
        xticklabels=True,
        yticklabels=True,
        #mask=mask_df.T,
        ax=ax)

        ax.set_title(f"{trimers_origin} Trimers KDLe Distance")
        ax.set_xlabel("Trimers")
        ax.set_ylabel("Repertoire")

        if hide_heamapx:
            ax.set_xticklabels([])
            ax.set_xticks([])

        if save_fig == True:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_figures\\LPA_heatmap_{}.png'.format(time)

            if os.path.exists("trimers_figures") is False:
                os.mkdir("trimers_figures")

            fig.savefig(save_path, bbox_inches='tight')

    # Performing PCA on the LPA distances results.      
    def pca(self, 
            scale:bool = True):
        
        """
        scale : bool -> scaling the data (standard scaler) before performing principle component analysis.
        """

        pca_input = self.sig_df.T
        pca_intput_index = self.sig_df.T.index

        if scale:
            pca_input = StandardScaler().fit_transform(self.sig_df.T)
        
        pca_obj = skPCA(n_components=0.95, random_state=42)
        pca_obj.fit(pca_input)
        pca_evar = pca_obj.explained_variance_ratio_
        pca_output = pca_obj.transform(pca_input)

        self.pca_df = pd.DataFrame(data = pca_output,
                                   index = pca_intput_index,
                                   columns = [f"PC{i+1}" for i in range(0,pca_output.shape[1])])
        
        self.pca_var = pd.DataFrame({"PC":self.pca_df.columns, "Variance":pca_evar})
        
        return self.pca_df, self.pca_var # pca results and explained variance by PC
    
    # Plotting the PCA results.
    def pca_scatter(self,
                    n_pcs:int = 3,
                    label_index:int = None,
                    save_fig:bool = True):
        """
        n_pcs : int -> number of principle components to present.
        label_index : int -> present the data with different colors based on unique values of the label when splitted by ".", unique index is the
                             location from which the unique values will be selected.
        save_fig : bool -> save the figure.
        """
        
        #
        x = np.array([1, 2, 3, 4])
        y = np.array([1, 0.94, 0.92, 0.90])
        pf = np.polyfit(x, y, 2)
        eq_pf = np.poly1d(pf)

        self.pca_df["label_index"] = "no_index"
        if isinstance(label_index, int):
            self.pca_df["label_index"] = self.pca_df.index.str.split(".").map(lambda X : X[label_index])
        
        pc_range = range(1,n_pcs+1)
        pc_ranges = [(i-1,j-1) for i in pc_range for j in pc_range]
        n_subplots = int(len(pc_ranges)**0.5)
        #print(f"n_pcs: {n_pcs} \npc_range: {pc_range} \npc_ranges: {pc_ranges} \nn_subplots: {n_subplots}")

        fig, axs = plt.subplots(n_subplots,n_subplots, figsize=(n_subplots*5,n_subplots*5))
        unique_labels = self.pca_df["label_index"].unique()

        for rc in pc_ranges:
            pc_x = f"PC{rc[0]+1}"
            pc_y = f"PC{rc[1]+1}"

            for unique_label in unique_labels:
                index_label_cond = (self.pca_df.label_index == unique_label)

                if n_pcs > 1:
                    temp_axs = axs[rc[0], rc[1]]
                else:
                    temp_axs = axs

                temp_axs.scatter(self.pca_df.loc[index_label_cond ,pc_x],
                                        self.pca_df.loc[index_label_cond ,pc_y],
                                        s = n_subplots*20,
                                        alpha = 0.6,
                                        edgecolor="black",
                                        lw = 1,
                                        label=unique_label)
                
                temp_axs.set_xlabel(pc_x, size=10)
                temp_axs.set_ylabel(pc_y, size=10)
                            
        ax_handles, ax_labels = temp_axs.get_legend_handles_labels() 
        fig.legend(handles=ax_handles,
                   labels=ax_labels,
                   loc='outside upper center', 
                   bbox_to_anchor=(0.5, eq_pf(n_pcs)),
                   fontsize=12,
                   ncol=len(unique_labels))
            
        if save_fig:
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            save_path = 'trimers_figures\\PCA_scatter_{}_{}.png'.format(label_index, time)

            if os.path.exists("trimers_figures") is False:
                    os.mkdir("trimers_figures")

            fig.savefig(save_path, bbox_inches='tight')

        plt.show()
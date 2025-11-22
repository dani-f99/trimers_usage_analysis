from datetime import datetime
import pandas as pd
import numpy as np
import os 
from tqdm import tqdm
import math

# Codons dictionary used for the nt sequence translation
codon_dic_updated = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': "S'", 'TCC': "S'", 'TCA': "S'", 'TCG': "S'",
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * for STOP
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S"', 'AGC': 'S"', 'AGA': 'R', 'AGG': 'R',

                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }

# Translating NT sequense and creating trimers dataframe
def nt_transalte_104(cdr_seq:str, 
                     aa_start:int = 1, 
                     aa_end:int = None):
                     #returns:str = "trimers",
                     #output_type:str = "pandas"):
    
    """
    cdr_seq: str -> cdr_seq.dna_seq format.
    aa_start: int -> first amino acid position.
    aa_end: int -> last aa_position.
    #returns: str -> non-functional argument.
    #output_type: str -> non functional argument.
    """
    
    # heavy chain varaible regions, see IMGT documentation for more information.
    regions_dict = {"fw1": np.arange(1,27),
                    "cdr1": np.arange(27,39),
                    "fw2": np.arange(39,56),
                    "cdr2": np.arange(56,66),
                    "fw3": np.arange(66,105)}
    
    # the first column should be the cdr3_aa and the second germline\seq
    cdr_seq = cdr_seq.split(".") #sequence of the cdr, changes between sequences.
    cdr3_len = len(cdr_seq[0]) #cdr3 length
    nt_seq = cdr_seq[1] #1-104 sequence

    translated = []
    t_length = int((len(nt_seq)-len(nt_seq)%3)/3)

    #if the number of NT spacers in the NT sequence isn't equal to 3 there is a change in 
    #Reading frame and need to cheek the sequence 
    n_spaces = nt_seq.count("-")
    if  n_spaces % 3 != 0:
        raise Exception("Number of NT spacers dosent divide by 3, cheek sequence") 

    # Translating the NT sequence
    for i in range(1,t_length):
        codon = nt_seq[i*3-3:i*3]
      
        if codon in list(codon_dic_updated.keys()):
            aa = codon_dic_updated[codon]  
        else:
            aa = "-"
         
        translated.append(aa)

    # If no end was decided (aa_end argument) the code will itirate over all the dna sequence
    if aa_end is None:
        aa_end = t_length

    results_aa = pd.Series(data=translated[aa_start-1:aa_end], index=range(aa_start,aa_end))
    results_aa_cleaned = results_aa[results_aa != "-"]
    
    # output dataframe creation
    r_dict = {"first_aa":[], "trimer":[], "ncount":1}
    j_count = 1
    j_positions = []
    for j in range(1, len(results_aa_cleaned)-1):
        r_dict["trimer"].append("".join(results_aa_cleaned.values[j-1:j+2]))

        # giving different numring scheme to the J region beacuse the CDR3 is highly variable in it's length
        first_aa = results_aa_cleaned.index[j-1]
        seq_up2_j = 104 + cdr3_len

        if first_aa > seq_up2_j:
            r_dict["first_aa"].append(f"{j_count}j")
            j_positions.append(f"{j_count}j")
            j_count += 1

        else:
            r_dict["first_aa"].append(first_aa)

    trimer_result = pd.DataFrame(r_dict)

    # determining the region of cdr3 and j based on the cdr3 length
    regions_dict["cdr3"] =  np.arange(105,105+cdr3_len)
    regions_dict["j"] = j_positions

    # assigning each trimer it's location
    reg_list = []
    for value in trimer_result.first_aa.values:
        for key in regions_dict.keys():
            if value in regions_dict[key]:
                reg_list.append(key)
    
    trimer_result.insert(loc=1, 
                         column="region", 
                         value= reg_list)

    return trimer_result #, reg_list # results_aa_cleaned

#
class Trimer():
    # list of the relevent columns of the sequences table sql output.
    relevent_cols = ["seq_id", "ai", "sample_id", "subject_id", "clone_id", "functional", "copy_number",  "cdr3_aa", "sequence", "germline"]
  
    # Initiating the trimer class
    def __init__(self, 
                 metadata_list : list,
                 seqs_loc = "trimers_raw_tables\\sequences.csv",
                 seq_collapsed_loc = "trimers_raw_tables\\sequence_collapse.csv",
                 metadata_loc = "trimers_raw_tables\\sample_metadata.csv",
                 rename_metadata = False,
                 new_metadata_names = None):
        
        '''
        metadata_list: array-like -> List of required metadata. Needed to join the relevent metadata from metadata SQL table to the sequcnes SQL table.
        seqs_loc: str -> Location of the sequences csv, can change it's location (defualt location is the import location of the ImportSQL class.
        metadata_loc: str -> Location of the metadata_table csv, can change it's location (defualt location is the import location of the ImportSQL class).
        rename_metadata: bool -> If rename of the metadata columns is requires, see new_metadata_names explanation).
        new_metadata_names: array-like -> Rename the column names of the import metadata (for example from "cell_subset_type" to "ab_target"). The list should contain
                                          the new names only in the same index rder as metadata_list input.
        '''

        # If the processed sequcnes csv already exists in the defualt location it will be imported.
        cleaned_seqs_path = "trimers_processed_tables\\cleaned_seqs.csv"
        processed_folder = "trimers_processed_tables"
        if os.path.exists("trimers_processed_tables\\cleaned_seqs.csv"):
            print("> Found cleaned_seqs.csv in the procesed tables folder.")
            self.cleaned_seqs = pd.read_csv(cleaned_seqs_path, index_col=0)
            print("> 'cleaned_seqs.csv' loaded.")

        # Creating the processed sequcnes file.
        else:
            # Importing the sequences csv
            # Filtring out: non-functional seqs, non-sample specific and  non-clone specific
            seqs = pd.read_csv(seqs_loc, index_col=0)[self.relevent_cols]
            self.seqs_col = pd.read_csv(seq_collapsed_loc, index_col=0)
            
            # Importing the metadata csv and orginizing the dataframe for the relevent information.
            metadata = pd.read_csv(metadata_loc, index_col=0)
            metadata_df = metadata.groupby(["sample_id","key"]).describe().reset_index()[[("sample_id",""), ("key",""), ("value","top")]].droplevel(level=1,axis=1)
            metadata_df = metadata_df[metadata_df.key.isin(metadata_list)]
            
            # Creating metadata dataframe in order to join it's values to the sequences dataframe
            data_mtdata = []
            for j in metadata_df.key.unique():
                data_mtdata.append(metadata_df[metadata_df.key==j].drop("key",axis=1).rename({"value":j},axis=1).reset_index(drop=True))

            result_metadata = pd.concat(data_mtdata, axis=1).T.drop_duplicates().T[["sample_id"] + metadata_list]

            # Renaming the metadata columns names according to the rename_metadata & new_metadata_names arguments
            if rename_metadata:
                rename_dict = {i:j for i,j in zip(metadata_list, new_metadata_names)}
                result_metadata.rename(rename_dict, axis=1, inplace=True)
                metadata_list = new_metadata_names
            else:
                metadata_list = metadata_list
            
            # Placing the metadata values into the sequcnes dataframe
            result_dict = {i[1]:[i[2],i[3]] for i in result_metadata.itertuples()}
            seqs[metadata_list] = list(seqs.sample_id.apply(lambda X : result_dict[X]).values)
            self.cleaned_seqs = seqs.copy()

            # Creating unique sequences only dataframe
            unique_seq_list = self.seqs_col.loc[(self.seqs_col.seq_ai.isin(self.cleaned_seqs.ai.values)) & (self.seqs_col.instances_in_subject > 0), "collapse_to_subject_seq_id"].values
            self.cleaned_seqs = self.cleaned_seqs[self.cleaned_seqs.seq_id.isin(unique_seq_list)]

            # Filtring out rows woth NT sequence that dosent divide by 3 and getting report df
            self.dropped_3dv = self.cleaned_seqs[(self.cleaned_seqs.germline.str.count("-")%3 != 0) | (self.cleaned_seqs.sequence.str.count("-")%3 != 0)]
            self.cleaned_seqs = self.cleaned_seqs[(self.cleaned_seqs.germline.str.count("-")%3 == 0) & (self.cleaned_seqs.sequence.str.count("-")%3 == 0)]

            # Saving the the cleaned sequences data into the defualt location
            if os.path.exists(processed_folder) == False:
                os.mkdir(processed_folder)
                print("> Processed folder havent been found, creating folder.")

            self.cleaned_seqs.to_csv(cleaned_seqs_path)
            ("> 'cleaned_seqs.csv' saved to processed_data folder.")
    
    # Creating trimmer dataframe for all the sub-datasets
    def create(self,
                      subdatasets_list:list,
                      source:str,
                      start_pos:int = None,
                      end_pos:int = None,
                      save_csv=True):
        """
        subdatasets_list: array-like -> list of the relevent columns in order to divide the dataset into sub-datasets.
        source: str -> Source of the DNA sequence:
                       a. 'germline' - will create trimmers based on the germline sequences of the unique clones.
                       b. 'top_seq' - will create trimmers based on the germline sequences of the top clone sequence (seq with the highers number of reads whitin clone).
                       c. 'all_seqs' - will create trimmers based on all of the sequences whitin a clone
        start_pos: int -> Starting aa position of the analysis.
        end_pos: int -> Final aa position of the anlysis.
        save_csv: bool -> Save the results? True by defualt.
        """
        
        # Creating labels for each sequence based on the subdatasets_list argument.
        labels = None
        for col in subdatasets_list:
            temp_val = self.cleaned_seqs[col].astype("str").str.replace(" ","")

            if labels is None:
                labels = temp_val
            else:
                labels += "." + temp_val

        # assiging the labels to the labels column
        self.cleaned_seqs["label"] = labels
       
        source_dict = {"germline":"germline", "top_seq":"sequence", "all_seq":"sequence"}
        if source in ["germline", "top_seq"]:
            idx = self.cleaned_seqs.groupby('clone_id')['copy_number'].idxmax()
            input_df = self.cleaned_seqs.loc[idx, ["label", "clone_id", "cdr3_aa"] + [source_dict[source]]]
        
        elif source == "all_seq":
            input_df = self.cleaned_seqs[["label","clone_id","cdr3_aa","sequence"]]

        else:
            raise Exception("Invalid source argument input string") 
        
        output_dfs = []
        unique_clones = []
        n_ulbl = len(self.cleaned_seqs.label.unique())
        print(f"> Itirating over sub-datasets (n={n_ulbl})")
        for ulbl in tqdm(self.cleaned_seqs.label.unique(), unit="dataset", initial=0):
            temp_output = input_df[input_df.label == ulbl]
            unique_clones.append([ulbl, temp_output.shape[0]]) # save the unique number of sequences per dataset
            temp_output = (temp_output.cdr3_aa + "." + temp_output[source_dict[source]]).apply(nt_transalte_104)
            temp_df = None

            #print(f"Concatenating {ulbl} sub-dataset analysis results") #[{}/{n_ulbl}]")
            for nij in temp_output:
                if nij is None:
                    temp_df = nij
                else:
                    temp_df = pd.concat([temp_df, nij], axis=0)
            
            temp_df["label"] = ulbl
            output_dfs.append(temp_df)

        self.trimers =pd.concat(output_dfs).groupby(["label","first_aa","region","trimer"]).count().reset_index()
        
        if save_csv:
            scsv_path = f"trimers_data\\{"-".join(subdatasets_list)}-{source}" #-trimers.csv"

            if os.path.exists("trimers_data") is False:
                os.mkdir("trimers_data")    
    
            time = datetime.now().strftime("[%d.%m.%y-%H;%M]")
            self.trimers.to_csv(scsv_path + f"-trimers_{time}.csv") #saves trimers output
            pd.DataFrame(unique_clones, columns=["label", "unique_seqs"]).to_csv(scsv_path + "-unqseqs.csv") #saves unique clones of datasets
            print("> Trimers and unique clones data saved to 'trimers_output' folder.")
            
        return self.trimers
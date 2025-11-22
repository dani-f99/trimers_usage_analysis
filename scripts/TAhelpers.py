import numpy as np
import pandas as pd
import math

# aa codon table
codon_dic_updated = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': "S", 'TCC': "S", 'TCA': "S", 'TCG': "S",
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * for STOP
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }

# Translating nt sequence
def nt_transalte_104(nt_seq):
    rem3 = len(nt_seq) % 3

    if rem3 != 0:
        nt_seq = nt_seq[:-rem3]

    translated = []
    seq_len = len(nt_seq)

    for i in range(1, int((seq_len+3)/3)):
        codon = nt_seq[i*3-3:i*3]
    
        if codon in list(codon_dic_updated.keys()):
            aa = codon_dic_updated[codon]  
        else:
            aa = "-"
            
        translated.append(aa)
    
    return translated



# custom function that match aa position to region
regions_dict = {"fw1": np.arange(1,27),
                "cdr1": np.arange(27,39),
                "fw2": np.arange(39,56),
                "cdr2": np.arange(56,66),
                "fw3": np.arange(66,105)}

# dictionary of regions positions
def find_region(pos_string:str) -> str:

    """
    pos_string : str -> amino acid position in string format
    """

    try:
        pos_int = int(pos_string)
        region = "cdr3"

        for r_key in regions_dict.keys():
            if pos_int in regions_dict[r_key]:
                region = r_key
                
    except:
        region = "j"
   
    return region


# Custom function for the processing of list of labels with the objective of getting unique values for each label index
def split_label_list(list_labels : list,
                     n_indexes : int) -> list:
    
    """
    this function is applied on a series of labels.
    list_labels : list -> list of labels to be splitted into it's values.
    n_indexes : int -> number of sub-labels in the joined label, divided by '.'
    """
    
    # If there are no labels returns nan value
    if (list_labels is np.nan):
        return [np.nan for i in range(0,n_indexes)]
    
    else:
        return [np.unique(np.array([i.split(".")[j] for i in list_labels])).tolist() for j in range(0,n_indexes)]
    

# Helper function that creates matrix of contidtions based on a dictionary, used for metadata filtring
def metadata_cond(dic_input: dict) -> pd.DataFrame:

    """
    metadata_cond : dict -> dictionary of label index:unique values (see line 162).
    """

    mt_dic = dic_input
    mt_keys = mt_dic.keys()
    mt_keys_len = np.array([len(mt_dic[i]) for i in mt_keys])

    n_columns = len(mt_keys)
    n_rows = np.prod(mt_keys_len)
    meta_df = pd.DataFrame(np.zeros((n_rows,n_columns)), columns = mt_keys)

    for i,v,l in zip(range(0,len(mt_keys)), mt_keys, mt_keys_len):
        array_length = n_rows
        unique_val = v
        
        if i == 0:          
            reps_numbers = np.prod(mt_keys_len[i+1:])

            temp_array = []
            for val_i in mt_dic[unique_val]:
                temp_array+=[val_i for i in range(1,reps_numbers+1)]

            meta_df[v] = temp_array

        else:
            reps_numbers = np.prod(mt_keys_len[i+1:])
            temp_array = []
            
            for val_i in mt_dic[unique_val]:
                temp_array += [val_i for k in range(1,reps_numbers+1)]

            final_array = temp_array.copy()
            for num in range(1,int(array_length/len(temp_array))):
                final_array += temp_array.copy()

            meta_df[v] = final_array  
            
    return meta_df

# Custom function for the calculation of amino acid diversity
def div_calc(array:np.array,
             order:int) -> int:
    
    """
    array : array-like -> input array of amino acid presence fractions, sum equal to 1.
    order : int -> order of the diversity calculation, 0=richness, >1 = diversity.
    """

    # Requirements cheek for functionality.
    cond_verification = {"array_type error":[all(isinstance(i, (float, int)) for i in array) is False, "String value in input array"],
                         "array_value error":[isinstance(array, (list, array, tuple)) is False, "array not in list,tupole or array format"],
                         "order_type error":[(isinstance(order, (int)) is False), "order input not integer"],
                         "order_value error":[order < 0, "need to be integer equal or bigger than zero"]}
    
    if len(array) == 1:
        return 1
    
    elif len(array) == 0:
        return 0
    
    else:
        # Error stop incase of invalid argument.
        if any([ib[0] for ib in cond_verification.values()]):
            cond_True = [i for i,j in zip(cond_verification.keys(), cond_verification.values()) if j[0] is True]
            raise Exception(f"Incorrect arguments input, see : {cond_True}. Fix according to documantation.")
            
        # Divesity Order 0 aka richness.
        if order == 0:
            div = len(array)

        # Divesity Order 1
        elif order == 1:
            Pi = [k * math.log(k, math.e) for k in array]
            div = math.exp(-sum(Pi))

        # Divesity order > 1
        else:
            Pi = [k**order for k in array]
            div = (sum(Pi))**(1/(1-order))
            
        return round(div)

# Calculate the fractions of seires value and return dictionary of its index:fraction
def dict_val_fraction(input_seires : pd.Series) -> dict:
    
    val_counts = input_seires.value_counts() / input_seires.value_counts().sum()
    val_dict = {i:round(float(j),3) for i,j in zip(val_counts.index, val_counts.values)}
    
    return val_dict

# Splitting label columns to seperate columns based on the label composition
def label_2cols(input_df : pd.DataFrame,
                index_labels : dict,
                mapper : dict = None) -> pd.DataFrame:
    """
    Split the label into columns of metadata accordinly.
    raw_trimers:pd.DataFrame -> raw trimers input.
    index_labels:dict -> Dictionary with the index and metadata names that construct the label.
    mapper:dict -> Dictionary with index and additional dictionary of mapper for example: {0:{"Non-SpikeB":"sn", "Spike+MemB":"sp"}}.
    """
    
    input = input_df
    idx_vals = index_labels

    for idx in idx_vals:
        i_key = idx
        i_value = idx_vals[i_key]

        input[i_value] = input.label.str.split(".").apply(lambda X : X[i_key])

        try:
            input[i_value] = input[i_value].map(mapper[idx])
        except:
            pass
        
    return input

# Simple function to return string according to number posistion compared to value
def return_2std(value : int, 
                std : int) -> str:
    """
    value : int -> The value we want to know it's position
    std : int -> Target value we comparing to.
    """

    if value > std:
        return "head"
    elif value < std:
        return "tail"
    else:
        return "middle"


# function that finds quantile of value in series
def find_qnt(value : int, 
             quantile : pd.Series) -> float:
    """
    value : int -> the value we want to see it's quantile location.
    quantile : series of quantile of the vector.
    """

    for i in range(0,len(quantile))[:-1]:
        i_low = quantile.iloc[int(i)]
        i_high = quantile.iloc[i+1]

        if (i_low <= value) & (value <= i_high):
            return float(quantile.index[i+1])
        
# Function that calcualte the diversity of each unique trimer in document. 
# applied to raw_trimers.groupby(["label","first_aa"]).apply().
# uses scripts.TAhelpers div_calc.
def raw2div(trimers_raw : pd.DataFrame,
            div : int = 0) -> pd.DataFrame:
    """ 
    order : int -> order of the diversity to calcualte.
    """

    raw_input =  trimers_raw[["trimer", "ncount"]].set_index("trimer",drop=True).sort_values("ncount", ascending=False) #re-orgnizing the data
    count_array = [float(i) for i in raw_input.ncount.values / np.sum(raw_input["ncount"].values)]  #array of trimers in amino acid position
    div = div_calc(count_array , order=div) # diversity calculation

    return raw_input[:div]

# custom function for region assigment, depends on the find_region function.
def assign_region(pos_array : list) -> list:
    """
    pos_array : list -> array of amino acid position in int format (preferbly).
    """
    output = []

    for i in pos_array:
        output.append(find_region(int(i)))

    return output

# custom functions that divide row of the trimers_count() format into multiple rows if more then 1 amino acid position is present in unq_pos.
def split_rows(input_df:pd.DataFrame,
               div:int = 0) -> pd.DataFrame:
    """
    Splitting rows with multiple position (in the TrimerAnalysis.TrimerPresence() format)
    input_df : pd.DataFrame -> TrimerPresence dataframe.
    div : int -> diversity of iterest
    """
    
    input_df = input_df[(input_df[f"div{div}_npos"].notnull()) & (input_df[f"div{div}_npos"] > 0)]
    original_columns = list(input_df.columns)
    modified_columns = original_columns[:2] + [f'div{div}_npos', f'div{div}_unqpos', f'div{div}_pcount',f'div{div}_ncount'] + original_columns[-3:]
 
    final_output = []

    for i in input_df[modified_columns].iterrows():
        row_data = i[1]
        n_positions = row_data[f'div{div}_npos']
        
        if n_positions == 1:
            i_vals = list(row_data.values)
            i_vals[2:6] = int(i_vals[2]), i_vals[3][0], list(i_vals[4].keys())[0] , round(list(i_vals[4].values())[0] * i_vals[5])
            final_output.append(i_vals)

        else:
            for j in row_data[f"div{div}_unqpos"]:
                j_vals = list(row_data.values)
                pos_frac = row_data.iloc[4][j]
                pos_count = round(pos_frac * row_data.iloc[5])

                j_vals[2:6] = int(j_vals[2]), j, j ,pos_count
                final_output.append(j_vals)

    final_output = pd.DataFrame(data=final_output, columns=modified_columns)
    final_output.drop(labels=f"div{div}_pcount", axis=1, inplace=True)
    final_output["region"] = final_output[f"div{div}_unqpos"].apply(find_region)
                
    return final_output
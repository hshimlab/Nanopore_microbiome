'''
2022.08.05
Yunseol Park
'''
import pandas as pd
import os
import statistics
import code

def iterate_pd(microbe_df, wimp_df, ref_list=None, genus=False, inconclusive=False):
    '''
    Function that separates the WIMP data into different host-microbiome interaction types.

    Args:
        microbe_df (DataFrame): data frame that contains the interaction type
        wimp_df (DataFrame): data frame that contains the WIMP data
        ref_list (list): list of all the microorganisms whose interaction types have already been identified
                         (also used as an indicator to check if we are looking at the data for the first time or not)
                         (default: None)
        genus (boolean): boolean that indicates whether or not to look at only the genus level
                         (default: False)
        inconclusive (boolean): boolean that indicates whether or not the search is for inconclusive interaction types
                                (default: False)
    Returns:
        same_dict (dict): nested dictionary with species name as key and its information as value
    '''
    same_dict = {}
    # Loop through each row of microbe_df
    for row in range(0, len(microbe_df)):
        df = microbe_df.iloc[row]
        name = df['Name']       # Species name
        if not inconclusive:
            # If we do not look at the genus level, then the name should be longer than 1 word
            if not genus and len(name.split()) == 1:
                continue
            # And if we do look at the genus level only, then the name should not be longer than 2 words
            if genus and len(name.split()) >= 2:
                continue
        # Formating name for the parenthesis
        if '(' in name:
            name = name.replace('(','\(').replace(')','\)')
        # Select WIMP data that starts with the species name
        s_df = wimp_df[wimp_df['species'].str.contains('^'+name, case=False)]
        # Convert into dictionary
        add = s_df.set_index('species').transpose().to_dict()
        # Select WIMP data - some data starts with 'unclassified' together with genus name
        second_df = wimp_df[wimp_df['species'].str.contains('unclassified ' + name, case=False)]
        add.update(second_df.set_index('species').transpose().to_dict())
        # Select WIMP data - some data have the genus names encased in brackets
        final_df = wimp_df[wimp_df['species'].str.contains('\[{}\] {}'.format(name.split()[0], ' '.join(name.split()[1:])), case=False)]
        add.update(final_df.set_index('species').transpose().to_dict())
        # If ref_list is given, then we only take the ones that are not in ref_list
        if ref_list:
            add_keys = list(add.keys())
            for key in add_keys:
                if key == '-' or key.lower() in ['homo sapiens', 'bacteria', 'root', 'cellular organisms']:
                    del add[key]
                elif key in ref_list:
                    del add[key]
        # Add the species information the the dictionary
        same_dict.update(add)
    return same_dict

def find_undefined(defined_list, wimp_df):
    '''
    Function that specifically finds undefined host-microbiome interaction type.

    Args:
        defined_list (list): list of species names that already have interaction types
        wimp_df (DataFrame): data frame that contains the WIMP data
    Returns:
        undefined_dict (dict): nested dictionary with species name as key and its information as value
    '''
    # Convert WIMP data frame into dictionary
    wimp_dict = wimp_df.set_index('species').transpose().to_dict()
    # Only add to the undefined if it is not in the defined_list
    undefined_dict = {}
    defined_list = [i.lower() for i in defined_list]
    for name, value in wimp_dict.items():
        if name.lower() not in defined_list:
            undefined_dict[name] = value
    return undefined_dict

def search_df(wimp_df, class_df):
    '''
    Function that finds all the beneficial, harmful, commensal, inconclusive, and undefined
    host-microbiome interaction types and counts them.

    Args:
        wimp_df (DataFrame): data frame that contains the WIMP data
        class_df (DataFrame): data frame that contains the classification data
    Returns:
        bene (DataFrame): data frame that contains the beneficial microbes
        harm (DataFrame): data frame that contains the harmful microbes
        comm (DataFrame): data frame that contains the commensal microbes
        inco (DataFrame): data frame that contains the inconclusive microbes
        undefined_dict ():
    '''
    # Split the classification data frame into different interaction types
    bene_df = class_df.loc[class_df['Interaction Type']=='beneficial']
    harm_df = class_df.loc[class_df['Interaction Type']=='harmful']
    comm_df = class_df.loc[class_df['Interaction Type']=='commensal']
    inco_df = class_df.loc[class_df['Interaction Type']=='inconclusive']
    # Get the WIMP data for each interaction type
    bene = iterate_pd(bene_df, wimp_df)
    harm = iterate_pd(harm_df, wimp_df)
    comm = iterate_pd(comm_df, wimp_df)
    # Get WIMP data for each interaction type again
    ## This time, we add ref_list and genus parameters
    ## We look at organisms at a higher level (genus level)
    ref_list = list(bene.keys()) + list(harm.keys()) + list(comm.keys())
    bene.update(iterate_pd(bene_df, wimp_df, ref_list, True))
    harm.update(iterate_pd(harm_df, wimp_df, ref_list, True))
    comm.update(iterate_pd(comm_df, wimp_df, ref_list, True))
    # Get WIMP data for inconclusive interaction type
    ref_list = list(bene.keys()) + list(harm.keys()) + list(comm.keys())
    inco = iterate_pd(inco_df, wimp_df, ref_list=ref_list, inconclusive=True)
    # Get WIMP data for undefined
    defined_list = list(bene.keys()) + list(harm.keys()) + list(comm.keys()) + list(inco.keys())
    undefined_dict = find_undefined(defined_list, wimp_df)
    return bene, harm, comm, inco, undefined_dict

def write_csv(write_path, write_name, dict_list, save_list):
    '''
    Function that writes WIMP information of data as CSV files according to their samples and
    host-microbiome interaction types.

    Args:
        write_path (str): directory name to save the results (sample name)
        write_name (str): file name to save the results (sample name + interaction type)
        dict_list (list): list that contains dictionaries with WIMP information
        save_list (list): list that contains strings indicating the interaction type
                          for each element in dict_list
    Returns:
        sum_dict (dict): dictionary that combines dict_list and save_list (WIMP data + interaction type)
    '''
    sum_dict = {}
    # Make directories if they don't already exist
    if not os.path.exists(write_path):
        os.makedirs(write_path)
    if not os.path.exists(write_path + '/' + write_name):
        os.makedirs(write_path + '/' + write_name)
    # Get the WIMP data and interaction type
    for wd_dict, save in zip(dict_list, save_list):
        # Write to CSV
        df = pd.DataFrame(wd_dict).transpose()
        df.to_csv('{}/{}/{}_{}.csv'.format(write_path, write_name, write_name, save))
        # Save WIMP data to the interaction type
        sum_dict[save] = df
    return sum_dict

def square_brackets(g):
    # Remove square brackets if any
    if '[' in g:
        g = g[1:-1]
    return g

def count_species(write_path, write_name, wimp_df):
    '''
    Function that separates the species into different levels (genus, species, strain, virus)
    and counts the number of organisms in each level. Here, genus includes all organisms that
    are of higher level than genus.

    Args:
        write_path (str): directory name to save the results (sample name)
        write_name (str): file name to save the results (sample name + interaction type)
        wimp_df (DataFrame): data frame that contains the WIMP data
    Returns:
        len(genus), len(species), len(strain), len(virus): lengths of each level
    '''
    # Get all the species names into a list
    item_name = list(wimp_df['species'])
    genus, species, strain, virus = set(), set(), set(), set()
    # Files to write down all the species names to
    g_file = open('{}/{}/{}_genus.txt'.format(write_path, write_name, write_name), 'w')
    s_file = open('{}/{}/{}_species.txt'.format(write_path, write_name, write_name), 'w')
    str_file = open('{}/{}/{}_strain.txt'.format(write_path, write_name, write_name), 'w')
    v_file = open('{}/{}/{}_virus.txt'.format(write_path, write_name, write_name), 'w')
    # Loop over each organism
    for name in item_name:
        if name == '-' or name.lower() in ['homo sapiens', 'bacteria', 'root']:
            continue
        # Viruses
        if 'phage' in name or 'virus' in name:
            virus.add(name.lower())
            v_file.write(name + '\n')
            continue
        if 'endosymbiont' in name:
            name = name.split('endosymbiont')[0].rstrip()
            if name == '':
                continue
        split_n = name.split()
        # Genus level
        if len(split_n) == 1:
            # Change name to make them the same with other organisms
            if name == 'Cardinium':
                name = 'Candidatus Cardinium'
            genus.add(name.lower())
            g_file.write(name + '\n')
        # Some genus level organisms are formatted as 'unclassified ___'
        elif 'unclassified' in name.lower():
            genus.add(' '.join(split_n[1:]).lower())
            g_file.write(' '.join(split_n[1:]) + '\n')
        # Some phylogenetic groups (higher level than genus level)
        elif 'group' in name.lower():
            # If there are two phylums in group
            if '/' in split_n[0]:
                genus.add(split_n[0].split('/')[0].lower())
                genus.add(split_n[0].split('/')[1].lower())
                g_file.write(split_n[0].split('/')[0] + '\n')
                g_file.write(split_n[0].split('/')[1] + '\n')
            else:
                genus.add(split_n[0].lower())
                g_file.write(split_n[0] + '\n')
            # If the species name is specified, then also add species
            if len(split_n) > 2:
                species.add(' '.join(split_n[:-1]).lower())
                s_file.write(' '.join(split_n[:-1]) + '\n')
        # Phylogenetic name with subdivision goes to genus (multiple)
        elif 'subdivision' in name.lower():
            genus.add(split_n[0].split('/')[0].lower())
            genus.add(split_n[0].split('/')[1].lower())
            g_file.write(split_n[0].split('/')[0] + '\n')
            g_file.write(split_n[0].split('/')[1] + '\n')
        # Same for 'complex'es
        elif 'complex' in name.lower():
            genus.add(split_n[0].lower())
            g_file.write(split_n[0] + '\n')
            # If there are more than one species, add both
            if '/' in split_n[1]:
                s_split = split_n[1].split('/')
                species.add('{} {}'.format(split_n[0], s_split[0]).lower())
                species.add('{} {}'.format(split_n[0], s_split[1]).lower())
                s_file.write('{} {}'.format(split_n[0], s_split[0]) + '\n')
                s_file.write('{} {}'.format(split_n[0], s_split[1]) + '\n')
            # Remove anything in parenthesis
            elif '(' in split_n[-1]:
                species.add(' '.join(split_n[:-2]).lower())
                s_file.write(' '.join(split_n[:-2]) + '\n')
            else:
                species.add(' '.join(split_n[:-1]).lower())
                s_file.write(' '.join(split_n[:-1]) + '\n')
        # The following elif statments are for exceptions
        elif name in ['Candidatus Kinetoplastibacterium blastocrithidii', 'Candidatus Nanosynbacter lyticus', 'Candidatus Thioglobus autotrophicus', 'Candidatus Planktophila dulcis', 'Candidatus Desulforudis audaxviator', 'Candidatus Tremblaya princeps', 'Candidatus Purcelliella pentastirinorum', 'Candidatus Phytoplasma asteris', 'Candidatus Planktophila vernalis', 'Candidatus Methylopumilus turicensis', 'Candidatus Promineofilum breve', 'Candidatus Fonsibacter ubiquis', 'Candidatus Kinetoplastibacterium sorsogonicusi', 'Candidatus Vampirococcus archaeovorus', 'Candidatus Gullanella endobia']:
            genus.add(' '.join(split_n[:2]).lower())
            species.add(name.lower())
            g_file.write(' '.join(split_n[:2]) + '\n')
            s_file.write(name + '\n')
        elif name in ['Candidatus Koribacter versatilis Ellin345' 'Candidatus Sulcia muelleri PSPU', 'Candidatus Amoebophilus asiaticus 5a2', 'Candidatus Atelocyanobacterium thalassa isolate ALOHA', 'Candidatus Carsonella ruddii PV', 'Candidatus Symbiobacter mobilis CR', 'Candidatus Solibacter usitatus Ellin6076', 'Candidatus Cloacimonas acidaminovorans str. Evry', 'Candidatus Protochlamydia amoebophila UWE25']:
            genus.add(' '.join(split_n[:2]).lower())
            species.add(' '.join(split_n[:3]).lower())
            strain.add(name.lower())
            g_file.write(' '.join(split_n[:2]) + '\n')
            s_file.write(' '.join(split_n[:3]) + '\n')
            str_file.write(name + '\n')
        elif name in ['Candidatus Erwinia haradaeae', 'Candidatus Filomicrobium marinum', 'Candidatus Thiodictyon syntrophicum', 'Candidatus Liberibacter asiaticus']:
            genus.add(split_n[1].lower())
            species.add(name.lower())
            g_file.write(split_n[1] + '\n')
            s_file.write(name + '\n')
        elif name == 'Candidatus Mycoplasma haemolamae str. Purdue':
            genus.add(split_n[1].lower())
            species.add(' '.join(split_n[:3]).lower())
            strain.add(name.lower())
            g_file.write(split_n[1] + '\n')
            s_file.write(' '.join(split_n[:3]) + '\n')
            str_file.write(name + '\n')
        elif name in ['Candidatus Arthromitus']:
            genus.add(name.lower())
            g_file.write(name + '\n')
        elif name in ['Candidatus Arthromitus sp. SFB-rat-Yit']:
            genus.add(' '.join(split_n[:2]).lower())
            strain.add(name.lower())
            g_file.write(' '.join(split_n[:2]) + '\n')
            str_file.write(name + '\n')
        elif name == 'Escherichia coli O43 str. RM10042':
            genus.add(split_n[0].lower())
            species.add(' '.join(split_n[:2]).lower())
            strain.add(name.lower())
            g_file.write(split_n[0] + '\n')
            s_file.write(' '.join(split_n[:2]) + '\n')
            str_file.write(name + '\n')
        elif 'Buchnera aphidicola' in name:
            genus.add(split_n[0].lower())
            species.add(' '.join(split_n[:2]).lower())
            strain.add(name.lower())
            g_file.write(split_n[0] + '\n')
            s_file.write(' '.join(split_n[:2]) + '\n')
            str_file.write(name + '\n')
        # Biotype or serovar or serogroups
        elif 'biotype' in name or 'serovar' in name or 'serogroup' in name:
            genus.add(square_brackets(split_n[0].lower()))
            species.add(' '.join(split_n[:2]).lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(' '.join(split_n[:2]) + '\n')
        # Subspecies
        elif 'subsp.' in name:
            subs_split = name.split('subsp.')
            genus.add(square_brackets(split_n[0].lower()))
            species.add(subs_split[0].rstrip().lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(subs_split[0].rstrip() + '\n')
            # Strains if it exists
            if 'str.' in name:
                strain.add(name.lower())
                str_file.write(name + '\n')
            if len(subs_split[1].split()) > 1:
                strain.add(name.lower())
                str_file.write(name + '\n')
        # If species names are formatted as 'sp. ___'
        elif 'sp.' in name:
            genus.add(square_brackets(split_n[0].lower()))
            g_file.write(square_brackets(split_n[0]) + '\n')
            # Add strains (having sp. also indicate that it is a strain that is not identified)
            if 'str.' in name:
                species.add(name.split('str.')[0].rstrip().lower())
                strain.add(name.lower())
                s_file.write(name.split('str.')[0].rstrip() + '\n')
                str_file.write(name + '\n')
            else:
                strain.add(name.lower())
                str_file.write(name + '\n')
        # Strains
        elif 'str.' in name:
            genus.add(square_brackets(split_n[0].lower()))
            species.add(name.split('str.')[0].rstrip().lower())
            strain.add(name.lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(name.split('str.')[0].rstrip() + '\n')
            str_file.write(name + '\n')
        # Pathovars and biovar
        elif 'pv.' in name or 'bv.' in name:
            if 'pv.' in name:
                abbreviation = 'pv.'
            elif 'bv.' in name:
                abbreviation = 'bv.'
            genus.add(square_brackets(split_n[0].lower()))
            species.add(name.split(abbreviation)[0].rstrip().lower())
            strain.add(name.lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(name.split(abbreviation)[0].rstrip() + '\n')
            str_file.write(name + '\n')
        # For those with genus, species, and strains in name
        elif len(split_n) > 2:
            genus.add(square_brackets(split_n[0].lower()))
            species.add(' '.join(split_n[:2]).lower())
            strain.add(name.lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(' '.join(split_n[:2]) + '\n')
            str_file.write(name + '\n')
        # Those that have square brackets around the genus name
        else:
            genus.add(square_brackets(split_n[0].lower()))
            species.add(name.lower())
            g_file.write(square_brackets(split_n[0]) + '\n')
            s_file.write(name + '\n')
    g_file.close()
    s_file.close()
    str_file.close()
    return len(genus), len(species), len(strain), len(virus)

def make_summary(write_path, write_name, sum_dict, wimp_df):
    '''
    Function that calculates the summary statistics of WIMP data.

    Args:
        write_path (str): directory name to save the results (sample name)
        write_name (str): file name to save the results (sample name + interaction type)
        sum_dict (dict): dictionary that contains WIMP data as value and interaction type as key
        wimp_df (DataFrame): data frame that contains the WIMP data
    Returns:
        table (dict): dictionary that contains the summary statistics
    '''
    table = {}
    reads = []      # Number of reads
    seqlen = []     # Length of sequence
    # Count the number of genus, species, strain, and virus
    genus, species, strain, virus = count_species(write_path, write_name, wimp_df)
    # Add to table
    table['# genus'] = genus
    table['# species'] = species
    table['# strain'] = strain
    table['# virus'] = virus
    for key, value in sum_dict.items():
        # Add number of species of each interaction type
        table[key] = len(value)
        # Get number of read and sequence length
        if not value.empty:
            reads += list(value['reads'])
            seqlen += list(value['seqlen_mean'])
    # Caluclate the average
    table['reads'] = statistics.mean(reads)
    table['seqlen'] = statistics.mean(seqlen)
    return table

def write_summary(write_path, summary_table):
    sum_df = pd.DataFrame(summary_table).transpose()
    col = list(sum_df.columns)
    col.remove('reads')
    col.remove('seqlen')
    sum_df = sum_df.reindex(columns=col[:4]+['reads', 'seqlen']+col[4:])
    #code.interact(local = dict(globals(), **locals()))
    sum_df.to_csv(write_path + '/summary_all_samples.csv')

def find_species(wimp_path, class_path, sample=['Stool', 'Saliva']):
    '''
    Function that classifies the species into host-microbiome interaction types, samples,
    and that generates a summary of the interaction.

    Args:
        wimp_path (str): directory name that contains all the WIMP files
        class_path (str): directory name that contains all the classification files
        sample (list): a list that indicates which sample to focus on.
                       (default: ['Stool', 'Saliva'])
    Returns:
        None
    '''
    # Path for the output
    write_path = 'results/' + class_path.split('/')[-1]
    # Create empty dictionary to store the summary statistics
    summary_table = {}
    for s in sample:
        # Get all file paths for WIMP data
        wimp_data = [i for i in sorted(os.listdir(wimp_path)) if s in i]
        print(wimp_data)
        # Get all file path for the sample
        class_data = [i for i in sorted(os.listdir(class_path)) if s in i][0]
        print(class_data)
        # Convert the classification file into a data frame
        class_df = pd.read_csv(class_path + '/' + class_data)
        for wd in wimp_data:
            # Convert the WIMP file into a data frame and clean up
            wimp_df = pd.read_csv(wimp_path + '/' + wd)
            wimp_df.rename(columns={'Unnamed: 0':'species'}, inplace=True)
            # Get information on beneficial, harmful, commensal, inconclusive, undefined species
            bene, harm, comm, inco, undefined_dict = search_df(wimp_df, class_df)
            # Write the information to CSV files
            dict_list = [bene, harm, comm, inco, undefined_dict]
            save_list = ['beneficial', 'harmful', 'commensal', 'inconclusive', 'undefined']
            write_name = '_'.join(wd.split('.')[0].split('_')[1:])
            sum_dict = write_csv(write_path, write_name, dict_list, save_list)
            # Calculate summary statistics
            table = make_summary(write_path, write_name, sum_dict, wimp_df)
            summary_table[write_name] = table
    write_summary(write_path, summary_table)

if __name__ == '__main__':
    wimp_path = '/home/yunseol/microbiome/results/WIMP'
    class_path = '/home/yunseol/microbiome/data/Classification'
    find_species(wimp_path, class_path)
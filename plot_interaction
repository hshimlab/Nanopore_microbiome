'''
2023.04.16
Yunseol Park
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def read_csv(filename):
    '''
    Function that reads in the summary CSV file and generates a data frame.
    
    Args:
        filename (str): name of the summary CSV file
    Returns:
        df (DataFrame): data frame that contains the information of different host-microbiome interaction types
    '''
    # Read in the file and convert into DataFrame
    df = pd.read_csv(filename)
    # Remove unnecessary columns
    df.drop(columns=['# genus','# species','# strain','# virus','reads','seqlen','undefined'], inplace=True)
    # Format sample names and set to index
    df['Unnamed: 0'] = df['Unnamed: 0'].str.split('_').str.join('\n')
    df = df.set_index('Unnamed: 0')
    df.index.name = None
    # Convert into percentage
    df = df.div(df.sum(axis=1), axis=0) * 100
    return df

def plot_stacked(filename, output):
    '''
    Function that plots the host-microbiome interaction types.

    Args:
        filename (str): name of the summary CSV file
        output (str): name of the image file to save
    Returns:
        None
    '''
    palette = sns.color_palette('pastel')[:4]      # Seaborn color palette
    # Read in the CSV file and convert into data frame
    df = read_csv(filename)
    # Plot and set parameters
    plt.rcParams["figure.figsize"] = [8, 6]
    df.plot.bar(stacked=True, rot=0, color=palette)
    plt.title("Percentage of microbiome classification")
    plt.ylabel("Percentage (%)")
    plt.xlabel("Samples")
    plt.legend(loc='upper right')
    # Save figure
    plt.savefig(output)

if __name__ == '__main__':
    summary = 'results/Classification/summary_all_samples.csv'
    plot_stacked(summary, 'results/image/plot_percentage_col2.png')
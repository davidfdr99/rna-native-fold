import os
import argparse 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

FILE_DIR = './PseudoEnergy'
SAVE_DIR = './InteractionProfiles'

def create_interaction_plots(save_dir: str, file_dir = FILE_DIR):
    """
    Creates pairwise interaction profile plots in the specified 'save_dir'. The default is the 'InteractionProfiles' folder.
    The x-axis represents the pseudoenergy and the y-axis the distances from 0-20 Å.
    """

    # Color palette for plots using seaborn package
    pal = sns.color_palette("deep", 10)

    counter = 0
    # Create Interaction Profile for each plot
    for filename in os.listdir(file_dir):

        f = os.path.join(file_dir, filename)
        # Checking if it is a file
        if os.path.isfile(f) and filename.endswith(".txt"):
            with open(f, "r") as i:
                y_axis = []
                for line in i:
                    y_axis.append(float(line.split()[0]))
                x_axis = list(range(21))
                
                title = filename.strip(".txt")
                
                plt.plot(x_axis, y_axis, color=pal[counter])
                plt.title(title, fontsize=12)
                plt.xlabel('Distance in Å', fontsize=11)
                plt.ylabel('Pseudo-Energy', fontsize=11)
                plt.grid(True)
                plt.savefig(f'{save_dir}/{title}.png')
                plt.close()

                counter += 1
        
    print(f'Plots successfully saved in {save_dir}')

def main(save_dir: str):
    """
    Takes the directory to store the files in, default is InteractionProfiles, and then creates the plots.
    """ 
    
    if save_dir != SAVE_DIR and not os.path.isdir(save_dir):
        try:
            os.mkdir(save_dir)
        except OSError:
            raise ("Failed to create directory for interaction plots. Please parse valid path or use default!")
    create_interaction_plots(save_dir=save_dir)

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('plot_dir', nargs='?', type=str, default=SAVE_DIR, help='Takes directory to save interaction plots in!')
    args=parser.parse_args()

    main(save_dir=args.plot_dir)
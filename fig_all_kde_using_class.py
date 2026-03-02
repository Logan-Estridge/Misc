import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

class KDEAnalyzer:
    """A tool for generating Kernel Density Estimate (KDE) plots for CV data.

    This class encapsulates the KDE plotting configurations for different types of CV 
    analyses (distances, torsions, etc.) from PLUMED and provides a standardized plotting interface.

    Attributes:
        species (str): The specific identifier for the CV (e.g., 'd8').
        conversion_factor (float): Multiplier to adjust units (e.g., 10.0 for nm to Å).
        title (str): The main title string for the generated plot.
        x_label (str): The label for the horizontal axis.
        save_name (str): The prefix used for the output image file.
        filename (str): The derived path to the input data file.
    """

    def __init__(self, species: str, prefix: str, conversion_factor: float, title: str, x_label: str, save_name: str) -> None:
        """Initializes the attributes for the KDEAnalyzer with instances corresponding to a particular CV's metadata.

        Args: 
            species: The species of CV, e.g. d8 for distance 8.
            prefix: Base name of the colvar.dat file for the particular CV. Used to locate the colvar file.
            conversion_factor: Scaling factor applied to the raw data values.
            title: Title of the plot corresponding to the particular CV.
            x_label: X axis label of the plot corresponding to the particular CV (including units).
            save_name: Base name for the saved PNG file.
        """
        self.species = species
        self.conversion_factor = conversion_factor
        self.title = title
        self.x_label = x_label
        self.save_name = save_name
        self.filename = f"{prefix}_{species}.dat" # derivative of {prefix} and {species} thus not passed into __init__

    @classmethod
    def distance(cls, species: str) -> "KDEAnalyzer":
        """Factory method for distance analysis.

        Args:
            species: The distance CV's identifier (e.g. 'd1')

        Returns:
            KDEAnalyzer: An instance configured with distance-specific units and labels for plotting.
        """
        return cls(
            species=species,
            prefix="distances",
            conversion_factor=10,
            title="Distance",
            x_label="Distance (Å)",
            save_name="dist_kde"
        )

    @classmethod
    def torsion(cls, species: str) -> "KDEAnalyzer":
        """Factory method for torsion analysis.

        Args:
            species: The torsion CV's identifier (e.g. 't1')

        Returns:
            KDEAnalyzer: An instance configured with torsion-specific units and labels for plotting.
        """
        return cls(
            species=species,
            prefix="torsions",
            conversion_factor=1,
            title="Torsion Angle",
            x_label="Angle (Radians)",
            save_name="torsions_kde"
        )

    def run(self) -> None:
        """Executes the data processing and plotting pipeline.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        Loads the PLUMED data from the derived filename, applies unit conversions, 
        and generates a formatted KDE plot. Skips execution if the file is missing 
        or contains invalid data.
        """
        # 1. Check if file exists
        if not os.path.exists(self.filename):
            print(f"Warning: File {self.filename} not found. Skipping {self.species}.")
            return

        try:
            # 2. Load Data
            df = pd.read_csv(self.filename, sep='\s+', comment='#', header=None)
            data = df.iloc[:, 1:].values * self.conversion_factor
            num_cols = data.shape[1]

            # 3. Plotting
            plt.figure(figsize=(10, 6))

            # We can still loop to create custom labels for the legend
            colors = sns.color_palette("colorblind", num_cols)
            for i in range(num_cols):
                sns.kdeplot(data[:, i], color=colors[i], label=f'Mol {i+1}', alpha=0.9)

            # 4. Formatting
            plt.title(f"{self.title} ({self.species}): {num_cols} mols")
            plt.xlabel(self.x_label)
            plt.ylabel("Density")

            if num_cols <= 20:
                plt.legend()

            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()

            # 5. Save
            plt.savefig(f"{self.save_name}_{self.species}.png", dpi=300)
            plt.close()

        except pd.errors.EmptyDataError:
            print(f"Error: {self.filename} contains no data.")
        except pd.errors.ParserError:
            print(f"Error: {self.filename} is poorly formatted (check your delimiters).")
        except ValueError as ve:
            print(f"Data Error: {ve}")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
        except PermissionError:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
            print(f"Permission Error: Cannot read {self.filename}.")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        except Exception as e:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            print(f"Unexpected Error on {self.species}: {e}")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# --- Execution ---                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
dist_species = [f"d{x}" for x in range(8, 9, 1)]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
tors_species = []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# Create task objects using the class methods                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
tasks = [KDEAnalyzer.distance(s) for s in dist_species] + \                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        [KDEAnalyzer.torsion(s) for s in tors_species]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# Run them                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
def main():                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    for task in tasks:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
        task.run()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
if __name__ == "__main__":                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    main()  

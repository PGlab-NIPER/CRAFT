# CRAFT_cmd
CRAFT_cmd, a command-line variant of the CRAFT tool, is designed for cavity prediction and the estimation of associated physicochemical properties.
CRAFT 
(Cavity Recognition Assisted by Flow Transfer algorithm)
User Manual
This user manual provides step-by-step instructions for installing Python, installing the SciPy library, downloading and unzipping the CRAFT tool files, downloading PDB files into the PDB folder, and running the algorithm on the command prompt. The algorithm provides identified cavities in descending order of volume of a cavity, along with information about cavity residue, residue number, atom number, atom type, chain, and fourteen physicochemical properties.

a) Install Python:
Visit the official Python website at https://www.python.org/downloads/
Choose the appropriate version of Python for your operating system and click on the download link.
Run the installer once the download is complete.
During installation, make sure to select the option to add Python to your PATH environment variable. This allows you to run Python from the command prompt.

b) Pip install SciPy:
Open the command prompt on your computer.
Type "pip install scipy" and press enter.
The installation process will begin and may take a few minutes to complete.
Once the installation is finished, you can start using the SciPy library in your Python projects.

c) Download the CRAFT tool file (CRAFT.zip) from the webpage:
Visit the webpage Link where the CRAFT tool file is available for download.
Click on the download link to download the CRAFT.zip file to your computer.

d) Unzip the downloaded file CRAFT.zip:
Locate the CRAFT.zip file on your computer and right-click on it.
Choose the option to "Extract" or "Extract all" to unzip the file.
Choose a destination folder to extract the files to.

e) Download PDB files into a PDB folder present inside the downloaded CRAFT folder:
Open the extracted CRAFT folder.
Locate the PDB folder inside it.
Download the PDB files that you want to analyze and save them to this PDB folder.
f) Open the command prompt inside the CRAFT folder:
Open the extracted CRAFT folder.
Hold down the Shift key and right-click inside the folder.
Choose the option to "Open command window here" or "Open PowerShell window here".

g) Use the command "Python main.py":
Type ```bash python main.py``` in the command prompt and press enter.
The algorithm will begin and take user input to start the cavity scan inside the given protein.

h) Follow the user instructions available on the command prompt:
Once the algorithm starts running, it will prompt you for the name of the PDB file that you want to analyze.
Enter the name of the PDB file and press enter.
The CRAFT tool provides information about the identified cavities in the PDB file, along with information about the cavity residue, residue number, atom number, atom type, chain, and fourteen physicochemical properties of a cavity.

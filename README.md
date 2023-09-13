# CRAFT_cmd
CRAFT_cmd, a command-line version of the CRAFT tool that is designed for cavity prediction and the estimation of associated physicochemical properties. The a web-based tool is also available at http://14.139.57.41:8080/

# CRAFT (Cavity Recognition Assisted by Flow Transfer algorithm)

* User Manual

This user manual provides step-by-step instructions for installing Python, installing the SciPy library, downloading and unzipping the CRAFT tool files, downloading PDB files into the PDB folder, and running the algorithm on the command prompt. The algorithm provides identified cavities in descending order of volume of a cavity, along with information about cavity residue, residue number, atom number, atom type, chain, and fourteen physicochemical properties.

1. Install Python:

Visit the official Python website at https://www.python.org/downloads/
Choose the appropriate version of Python for your operating system and click on the download link.
Run the installer once the download is complete.
During installation, make sure to select the option to add Python to your PATH environment variable. This allows you to run Python from the command prompt.

2. Pip install SciPy:

Open the command prompt on your computer.
Type ```pip install scipy``` and press enter.

The installation process will begin and may take a few minutes to complete.
Once the installation is finished, you can start using the SciPy library in your Python projects.

3. Pip install freesasa:

To install FreeSASA, open your computer's command prompt.
Enter the command ```pip install freesasa``` and press enter.
 
The installation process will commence and might require a few minutes to finalize. Once installed, you can readily incorporate the FreeSASA library into your Python projects.

4. Download the CRAFT tool file (CRAFT_cmd.zip):
   
Download the CRAFT_cmd.zip file from curretn github repository ().

5. Unzip the downloaded file CRAFT_cmd.zip:

Locate the CRAFT_cmd.zip file on your computer and right-click on it.
Choose the option to "Extract" or "Extract all" to unzip the file.
Choose a destination folder to extract the files to.

6. Download PDB files into a PDB folder present inside the downloaded CRAFT folder:

Open the extracted CRAFT folder.
Locate the PDB folder inside it.
Download the PDB files that you want to analyze and save them to this PDB folder.

7. Open the command prompt inside the CRAFT folder:
   
Open the extracted CRAFT folder.
Hold down the Shift key and right-click inside the folder.
Choose the option to "Open command window here" or "Open PowerShell window here".

8. Use the command "Python main.py":

Type ```python main.py``` in the command prompt and press enter.
The algorithm will begin and take user input to start the cavity scan inside the given protein.

9. Follow the user instructions available on the command prompt:

Once the algorithm starts running, it will prompt you for the name of the PDB file that you want to analyze.
Enter the name of the PDB file and press enter.
The CRAFT tool provides information about the identified cavities in the PDB file, along with information about the cavity residue, residue number, atom number, atom type, chain, and fourteen physicochemical properties of a cavity.

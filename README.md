# RenShape

<img src="https://github.com/matteo-brg/RenShape/blob/master/RenShape_logo.png">

**RenShape** (Reactor (anti)Neutrino spectral Shape) is a tool for evaluating the antineutrino spectrum from reactors using the summation method.

The (anti) neutrino spectra are evaluated using BetaShape (http://www.lnhb.fr/home/rd-activities/spectrum-processing-software/).

You can choose to use the database already compiled ("verza.lazy", October 2024), or compile your own.
For the latter, please refer to the README.md in the script folder.

For more information about the software, please take a look at the RenShape tutorial in the notebook folder.

### Installation

1. Open the terminal and clone the repository locally:
   Clone using a password protected SSH key (preferred)

   ```
   $ git clone 
   ```

   or Clone with HTTPS using the web URL

   ```
   $ git clone https://github.com/matteo-brg/RenShape
   ```
2. Install virtualenv using `python >= 3.8`

   ```
   $ sudo pip install virtualenv
   ```

   Create the lazyspectra virtual environment based on python3.X (X>=8)

   ```
   $ virtualenv -p `which python3.X` ~/lazyenv
   ```

   Run the following command to activate the virtual environment

   ```
   $ source lazyenv/bin/activate
   ```

   (optional) Check the version of python being used in the virtual environment

   ```
   (lazyenv) $ python -V
   ```
3. Install the required Python packages and dependencies. With the virtual environment active, move into the `RenShape` folder

   ```
   (lazyenv) $ cd RenShape
   ```

   and type

   ```
   (lazyenv) $ pip install .
   ```
4. (Optional:) Install Betashape following the instruction http://www.lnhb.fr/rd-activities/spectrum-processing-software/ (V2.4 tested).

   This step is required only if you want to create your own dataset or experiment with BetaShape itself (see the BetaShape tutorial inside notebook folder).

   I suggest you to follow these steps: download the betashape folder, then make every file inside the main betashape folder executable:

   ```
   chmod +x file_name
   ```

   I highly suggest you to create a virtual link to the betashape directory. If the betashape directory is located at "betashape_path", create a virtual link to it

   e.g.

   ```
   ln -s betashape_path betashape_rightversion
   ```

   This way, if you download a new version of betashape, you only need to update the link.

   Ensure that the BetaShape directory is defined in your PATH. Add the following lines in your .bashrc file

   ```
   #ADDED FOR BETASHAPE
   export PATH="$PATH:path_target/betashape_rightversion"
   ```

   reboot your terminal, and you should be good to go!

### Project top-level directory layout

    RenShape
    │
    ├── src                            # Project source code
    ├── scripts                        # Directory for scripts and executables
    ├── notebook                       # Directory for tutorials
    ├── requirements.txt               # Lists of packages to install
    ├── data                           # Directory for data
    ├── setup.py                       #
    ├── .gitignore                     # Specifies intentionally untracked files to ignore
    └── README.md                      # README file

 ASCII art tree structure taken from [here](https://codepen.io/patrickhlauke/pen/azbYWZ)

### Documentation

http://www.lnhb.fr/home/rd-activities/spectrum-processing-software/ BetaShape and its manual.

https://www.nndc.bnl.gov/ensdf/ ENSDF database

https://www.oecd-nea.org/dbdata/jeff/jeff33/ JEFF 3.3 database

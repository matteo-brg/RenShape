# HOWTO

### Unpack ENSDF file

2. Download the ENSDF files from a proper database (e.g. https://www.nndc.bnl.gov/ensdf/)

   Select the decay (B-) and the A range (e.g. 6-182)

   Press "select all" and "ENSDF text format".

   Save the resulting web page as a .txt (e.g. Selected_ENSDF_Datasets.txt).

   Unpack the txt file

   ```
   $ unpackENSDF.py -f data/SelectedENSDFDataset -o data/ENSDF
   ```

### Create a .lazy file

1. Create a database for RenShape. BetaShape must be installed in your system.

   Download the latest cumulative fission yields from the JEFF database (https://www.oecd-nea.org/dbdata/jeff/jeff33/ , Neutron Induced Fission Product Yields).

   Download the ENSDF data and unpack it (see previous script).

   ```
   $ createLazy.py -lp data/cavolo.lazy -jp data/JEFF33-nfy.txt -esp data/ENSDF -ep data/endfb-bm.txt -bp betashape_path -fix 1 
   ```
   **parameters:**

   -lp : string. Path where the .lazy file will be created. example: -lp data/cavolo.lazy. Default is None

   -jp : string. Path to the jeff database. example: data/JEFF33-nfy.txt Default is None.

   -jr : string. Optional. Version for Jeff. example: -jr 3.3. Defailt is "".

   -br : string. Optional. Version for BetaShape. example: -br 2.4. Default is "".

   -esp : string. Path to the ENSDF folder. Inside it, multiple .ensdf files are expected. example: data/ENSDF. Default is None.

   -ep : string. Path to the ENSDF/B sub-library. example: data/endfb-bm.txt. Default is None.

   -bp : string. Path to BetaShape directory. example: -bp home/BetaShape/betashape_rightversion. Default is None.

   -bo : string. Path for the temporary BetaShape output (where the dummyfolder101 will be created). See BetaShape tutorial for more details. If None, the current location will be used. Default is None.

   -bc : string. Option for BetaShape. See BetaShape manual for more details. "myEstep" and "nu" are required. Default is "myEstep=1 nu =1".

   -fix : int. Use ENSDF database to fix missing/continuum data. Default is 1.

   -ovr : int. Overwrite existing data. Default is 1.

   -t: str. Select "data_beta" for evaluating the electron spectra or "data_nu" for evaluating the neutrino spectra. Default is "data_nu"

# HOWTO

### Create a .lazy file
1. Create an hdf5 file (example.lazy) from the iaea Live Chart of Nuclides (https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html)
   ```
   $ extractNuchart.py -p /location/of_the/output_file -n outputfile_name -c config_file 
   ```
   ```
   $ extractNuchart.py -p /home/user/Desktop -n dummy -c /home/user/Desktop/lazyspectra/config/nuchart.cfg -not "This is an example" 
   ```

   **default parameters:** --note ('')
   it will create also a spectra.log file. Check it out!

2. Download the ENSDF files from a proper database (e.g. https://www.nndc.bnl.gov/ensdf/)

   Select the decay (B-) and the A range (e.g. 6-182)

   Press "select all" and "ENSDF text format".

   Save the resulting web page as a txt. (e.g. Selected_ENSDF_Datasets.txt). 

   Unpack the txt file
   ```
   $ unpackENSDF.py -f ENSDF_txt_filename -o output/directory/where_the_data_will_be_saved 
   ```   

3. Run BetaShape

   Change the betashape config file. In particular, change the path_to_betashape to your betashape directory

   Execute Betashape. All the files in the specified directory will be processed
   ```
   $ runBetashape.py -f .lazy_file -p path/of/the_ENSDF/directory -c config_file_to_use -o output/directory/for/betashape/output 
   ```   
   ```
   $ runBetashape.py -p /home/user/Desktop/dummy.lazy -p /home/user/ENSDF_dir -c /home/user/Desktop/lazyspectra/config/betashape.cfg -o /home/user/Desktop/betashape_output
   ```

4. Merge betashape output with the .lazy file
   ```
   $ mergeBetashapeLazy.py -f .lazy_file -p path/of/betashape/directory -t full
   ```   
   ```
   $ mergeBetashapeLazy.py -p /home/user/Desktop/dummy.lazy -p /home/user/Desktop/betashape_output -t full
   ```
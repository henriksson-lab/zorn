
# 
# Automatic:
#   
#   1. Run the command "download-db.sh" to automatically download and extract to:
#   /opt/software/conda_env/share/gtdbtk-2.4.0/db/
#   
#   Manual:
#   
#   1. Manually download the latest reference data:
#   wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
# 
# 2. Extract the archive to a target directory:
#   tar -xvzf gtdbtk_r220_data.tar.gz -C "/path/to/target/db" --strip 1 > /dev/null
# rm gtdbtk_r220_data.tar.gz
# 
# 3. Set the GTDBTK_DATA_PATH environment variable by running:
#   conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"
# 
# 
# - 
#   
#   
#   
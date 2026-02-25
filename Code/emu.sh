conda activate nanopack
conda activate emu

# EMU Wrapper available at https://github.com/marcelo-nd/EmuWrapper
EMUWRAPPER_LOC="/EmuWrapper"

### Screening

# Unzip sequences
export SC100_Scree=/1_Screening

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_Scree/Sequences -o $SC100_Scree

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_Scree/fastq -o $SC100_Scree -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_Scree/fastq_qc -o $SC100_Scree -d $EMU_DATABASE_DIR -c "TRUE" -p /databases/16s_copies.csv

### Selected SynComs Timepoints

# Unzip sequences
export SC100_TP=/2_Timepoints

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_TP -o $SC100_TP

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_TP/fastq -o $SC100_TP -q 10 -l 1000 -h 5000

# emu

export EMU_DATABASE_DIR=/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_TP/fastq_qc -o $SC100_TP -d $EMU_DATABASE_DIR -c "TRUE" -p /databases/16s_copies.csv


### Select SynComs

# Unzip sequences
export SC100_sel=/3_selected_repetition

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_sel -o $SC100_sel

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_sel/fastq -o $SC100_sel -q 10 -l 1000 -h 5000

# emu

export EMU_DATABASE_DIR=/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_sel/fastq_qc -o $SC100_sel -d $EMU_DATABASE_DIR -c "TRUE" -p /databases/16s_copies.csv


### Cocultures
# Unzip sequences
export SC100_sel=/4_Cocultures

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_sel -o $SC100_sel

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_sel/fastq -o $SC100_sel -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_sel/fastq_qc -o $SC100_sel -d $EMU_DATABASE_DIR -c "TRUE" -p /databases/16s_copies.csv


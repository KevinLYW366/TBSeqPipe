##### Common Module #####

#######################
# Load python library #
#######################
from snakemake.utils import min_version
import pandas as pd

###########################
# Check snakemake version #
###########################
min_version("6.13.1")

####################
# Read config file #
####################
configfile: "config/config.yaml"

################
# Load options #
################
# load sample list from sample_list.txt
with open(config["sample_list"], 'r') as f:
    SAMPLES = f.read().strip().split('\n')

# load fastq file name options
if config["fastq_read_id_format"]:
    FASTQREAD = "R"
else:
    FASTQREAD = ""

if config["fastq_suffix_format"]:
    FASTQSUFFIX = "fq.gz"
else:
    FASTQSUFFIX = "fastq.gz"

# load data directory format option
if config["data_dir_format"]:
    DATADIRINDEX = "/{sample}"
else:
    DATADIRINDEX = ""

##### Helper functions #####
def extract_dirname(wildcard, file_path):
    """
    Extract directory path from a file path.
    Example usage:
    Directory path extracted from input or output section could be used in params section.
    """
    dir_name = os.path.dirname(file_path)
    return dir_name


def get_fastq_reads(sample):
    """
    If data cleaning step is chosen to skip, just use the cleaned fastq files as the input of bwa_map
    """
    if config["data_cleaning"]:
        return ["results/00.data_cleaning/%s/%s_1.clean.fastq.gz" % (sample, sample),
                "results/00.data_cleaning/%s/%s_2.clean.fastq.gz" % (sample, sample)]
    else:
        if config["data_dir_format"]:
            return [config["data_dir"] + "/%s/%s_%s1.%s" % (sample, sample, FASTQREAD, FASTQSUFFIX),
                    config["data_dir"] + "/%s/%s_%s2.%s" % (sample, sample, FASTQREAD, FASTQSUFFIX)]
        else:
            return [config["data_dir"] + "/%s_%s1.%s" % (sample, FASTQREAD, FASTQSUFFIX),
                    config["data_dir"] + "/%s_%s2.%s" % (sample, FASTQREAD, FASTQSUFFIX)]


def get_qualified_results(path, wildcards):
    """
    Input a file path with sample wildcards such as "results/01.Alignment/{sample}/{sample}.sorted.bam".
    Return required file path of qualified samples based on kraken qc result.
    This function will be used to generate an output file list, which is helpful to
      prevent unqualified samples from running specific analysis steps of the workflow.
    """
    thres = config["kraken_threshold"]
    if config["kraken_filter"]:
        qc = pd.read_csv(checkpoints.kraken_qc.get().output[0], sep="\t", header=None, comment='#')
        qc.columns = ["sample", "percent", "category"]
        return expand(path, sample=qc[qc["percent"] > thres]["sample"])
    else:
        return expand(path, sample=SAMPLES)


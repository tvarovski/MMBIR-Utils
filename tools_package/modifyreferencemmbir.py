import regex as re
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s:.%(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

PATH="Mus_musculus.GRCm39.dna.primary_assembly.fa"
OUTPUT_FILE_NAME="out.fa"

def extractSeqFromFileToFileMM(file_path_in: str, file_path_out: str) -> None:
    """
    Extracts the DNA sequence from a reference genome file and saves it to a new file in a format that can be used by the MMBIR pipeline.
    Good for Mouse Genome Reference.

    Args:
        file_path_in (str): The path to the input reference genome file.
        file_path_out (str): The path to the output file where the DNA sequence will be saved.

    Returns:
        None: This function does not return anything, it only saves the DNA sequence to a file.

    Example:
        extractSeqFromFileToFileMM("reference_genome.fasta", "chromosome1.fasta")
    """
    #good for mouse genome
    file_in = open(file_path_in, "r")

    save=False
    readlines=True
    linecounter=0
    while readlines:
        linecounter+=1
        try:
            line = file_in.readline()
        except:
            logging.error(f"Failed to read line {linecounter}. EOF? Exiting")
            readlines=False
            break
        
        try:
            if line[0] == ">":
                if re.search(r"^>\d+|MT|[XY] dna:chromosome chromosome:GRCm39:[1-9XYM]+:\d+:\d+:\d+ REF$", line) != None:
                    logging.info(f"Extracting at line {linecounter}: {line.strip()}")
                    linewords = line.split()
                    chromosome=linewords[0].strip(">").strip()
                    if (chromosome=="1") | (chromosome=="2"):
                        chromosome=f">chr0{chromosome}\n"
                    elif chromosome=="X":
                        chromosome=f">chrX\n"
                    elif chromosome=="Y":
                        chromosome=f">chrY\n"
                    elif chromosome=="MT":
                        chromosome=f">chrM\n"
                    else:
                        chromosome=f">chr{chromosome}\n"
                    line=chromosome
                    logging.info(f"Found {line}")
                    save=True
                #fix for mitochondrial chr
                #elif re.search("mitochondrion, complete genome$", line) != None:
                #    print(line)
                #    line=f">chrM\n"
                #    print(line)
                #    save=True
                else:
                    save=False
        except:
            logging.error(f"Failed to parse line {linecounter}: {line}")
            raise Exception(f"Failed to parse line")
            


        if save:
            with open(file_path_out, "a") as file_out:
                file_out.write(line)

def extractSeqFromFileToFileHG(file_path_in: str, file_path_out: str) -> None:
    """
    Extracts the DNA sequence from a reference genome file and saves it to a new file in a format that can be used by the MMBIR pipeline.
    Good for Human Genome Reference.

    Args:
        file_path_in (str): The path to the input reference genome file.
        file_path_out (str): The path to the output file where the DNA sequence will be saved.

    Returns:
        None: This function does not return anything, it only saves the DNA sequence to a file.

    Example:
        extractSeqFromFileToFileMM("reference_genome.fasta", "chromosome1.fasta")
    """
    file_in = open(file_path_in, "r")

    save=False
    readlines=True
    while readlines:
        try:
            line = file_in.readline()
        except:
            logging.error("Failed to read line. EOF? Exiting")
            readlines=False
            break

        if line[0] == ">":
            if re.search("chromosome.*Primary.Assembly$", line) != None:
                logging.info(f"Extracting from: {line}")
                linewords = line.split()
                chromosome=linewords[4].strip(",")
                if (chromosome=="1") | (chromosome=="2"):
                    chromosome=f">chr0{chromosome}\n"
                elif chromosome=="X":
                    chromosome=f">chrX\n"
                elif chromosome=="Y":
                    chromosome=f">chrY\n"
                else:
                    chromosome=f">chr{chromosome}\n"
                line=chromosome
                logging.info(f"Found {line}")
                save=True
            #fix for mitochondrial chr
            elif re.search("mitochondrion, complete genome$", line) != None:
                logging.debug(line)
                line=f">chrM\n"
                logging.debug(line)
                save=True
            else:
                save=False

        if save:
            with open(file_path_out, "a") as file_out:
                file_out.write(line)

if __name__ == "__main__":

    extractSeqFromFileToFileMM(PATH, OUTPUT_FILE_NAME)
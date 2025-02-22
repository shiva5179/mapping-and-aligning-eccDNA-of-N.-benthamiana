import subprocess
import argparse
import os
import logging

# Configure logging for better verbosity and debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(cmd, log_prefix):
    """
    Executes a shell command and logs its output with detailed error handling.

    Args:
        cmd (list): Command to execute as a list of strings.
        log_prefix (str): Prefix for log messages to identify the command source.

    Returns:
        int: Return code of the command (0 for success, non-zero for failure).
    """
    cmd_str = ' '.join(cmd)
    logging.info(f"{log_prefix}: Running command: {cmd_str}")

    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   universal_newlines=True)
        stdout, stderr = process.communicate()

        if stdout:
            logging.info(f"{log_prefix}: Standard Output:\n{stdout}")
        if stderr:
            logging.error(f"{log_prefix}: Standard Error:\n{stderr}")

        return process.returncode

    except FileNotFoundError as e:
        logging.error(f"{log_prefix}: Command not found: {e}. Ensure it is in your PATH.")
        return 1  # Command not found is a failure

    except Exception as e:
        logging.error(f"{log_prefix}: Error executing command: {e}")
        return 1  # Generic error


def unmap_reads(input_file, output_file, reference_genome, log_prefix):
    """
    Maps reads to a reference genome using Minimap2, filters unmapped reads using Samtools,
    and sorts the output.

    Args:
        input_file (str): Path to the input FASTQ/BAM file.
        output_file (str): Path to the output BAM file containing unmapped reads.
        reference_genome (str): Path to the reference genome FASTA file.
        log_prefix (str): Prefix for log messages.

    Returns:
        bool: True if the process succeeds, False otherwise.
    """

    logging.info(f"{log_prefix}: Starting unmapping process...")

    # Construct commands
    minimap2_command = [
        "minimap2",
        "-ax", "sr",
        "-t", str(os.cpu_count()),
        reference_genome,
        input_file
    ]

    samtools_view_command = [
        "samtools",
        "view",
        "-b",
        "-f", "4",  # Changed from 12 to 4 to filter for unmapped reads
        "-"
    ]

    samtools_sort_command = [
        "samtools",
        "sort",
        "-o", output_file
    ]

    # Execute pipeline
    try:
        minimap2_process = subprocess.Popen(minimap2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        samtools_view_process = subprocess.Popen(samtools_view_command, stdin=minimap2_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        samtools_sort_process = subprocess.Popen(samtools_sort_command, stdin=samtools_view_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        # Close standard input for subprocesses to avoid deadlocks
        minimap2_process.stdout.close()
        samtools_view_process.stdout.close()

        # Capture output/error logs
        stdout, stderr = samtools_sort_process.communicate()

        if stdout:
            logging.info(f"{log_prefix}: Standard Output:\n{stdout}")
        if stderr:
            logging.error(f"{log_prefix}: Standard Error:\n{stderr}")

        # Check return codes
        if minimap2_process.returncode != 0 or samtools_view_process.returncode != 0 or samtools_sort_process.returncode != 0:
            logging.error(f"{log_prefix}: One or more processes failed.")
            return False

        logging.info(f"{log_prefix}: Unmapping process completed successfully.")
        return True

    except Exception as e:
        logging.error(f"{log_prefix}: An error occurred during the unmapping process: {e}")
        return False


def convert_bam_to_fastq(input_bam, output_fastq, log_prefix):
    """
    Converts a BAM file to FASTQ format using samtools.

    Args:
        input_bam (str): Path to the input BAM file.
        output_fastq (str): Path to the output FASTQ file.
        log_prefix (str): Prefix for log messages.

    Returns:
        bool: True if the conversion succeeds, False otherwise.
    """

    logging.info(f"{log_prefix}: Converting BAM to FASTQ format...")

    bam_to_fastq_command = [
        "samtools",
        "fastq",
        "-o", output_fastq,
        input_bam
    ]

    return_code = run_command(bam_to_fastq_command, log_prefix)
    if return_code == 0:
        logging.info(f"{log_prefix}: Conversion to FASTQ completed successfully. Output saved to {output_fastq}.")
        return True
    else:
        logging.error(f"{log_prefix}: Error: Failed to convert BAM to FASTQ.")
        return False


def cleanup_intermediate_files(files_to_remove):
    """
    Removes intermediate files.

    Args:
        files_to_remove (list): List of file paths to remove.
    """

    logging.info("Cleaning up intermediate files...")
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
            logging.info(f"Removed {file}")
        else:
            logging.warning(f"File {file} not found, skipping removal.")


def main():
    """
    Main function to handle arguments and process the unmapped reads step-by-step.
    """

    # Parsing command-line arguments
    parser = argparse.ArgumentParser(description="Unmap reads from raw FASTQ file using Minimap2 and Samtools.")
    parser.add_argument("--input", required=True, help="Path to the input FASTQ file containing raw reads.")
    parser.add_argument("--output", required=True, help="Path to the output FASTQ file after unmapping.")
    parser.add_argument("--mito", required=True, help="Path to the mitochondrial reference genome FASTA file.")
    parser.add_argument("--chloro", required=True, help="Path to the chloroplast reference genome FASTA file.")
    parser.add_argument("--viral_adna", required=True, help="Path to the viral ADNA reference genome FASTA file.")
    parser.add_argument("--viral_bsat", required=True, help="Path to the viral beta-sat reference genome FASTA file.")

    args = parser.parse_args()

    # Define intermediate file names
    mito_unmapped_bam = "mito_unmapped.bam"
    chloro_unmapped_bam = "chloro_unmapped.bam"
    viral_adna_unmapped_bam = "viral_adna_unmapped.bam"
    final_unmapped_bam = "final_unmapped.bam"

    # Unmap reads step-by-step through different reference genomes
    if not unmap_reads(args.input, mito_unmapped_bam, args.mito, "Mito_Unmapping"):
        logging.error("Mitochondrial unmapping failed. Exiting.")
        return

    if not unmap_reads(mito_unmapped_bam, chloro_unmapped_bam, args.chloro, "Chloro_Unmapping"):
        logging.error("Chloroplast unmapping failed. Exiting.")
        return

    if not unmap_reads(chloro_unmapped_bam, viral_adna_unmapped_bam, args.viral_adna, "Viral_ADNA_Unmapping"):
        logging.error("Viral ADNA unmapping failed. Exiting.")
        return

    if not unmap_reads(viral_adna_unmapped_bam, final_unmapped_bam, args.viral_bsat, "Viral_BetaSat_Unmapping"):
        logging.error("Viral BetaSat unmapping failed. Exiting.")
        return

    # Convert final BAM file to FASTQ format
    if not convert_bam_to_fastq(final_unmapped_bam, args.output, "BAM_to_FASTQ"):
        logging.error("Conversion to FASTQ failed. Exiting.")
        return

    # Clean up intermediate BAM files
    cleanup_intermediate_files([mito_unmapped_bam, chloro_unmapped_bam, viral_adna_unmapped_bam, final_unmapped_bam])

    logging.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()

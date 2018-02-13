import argparse
import base64
import json
import logging as log
import os
import random
import shutil
from io import StringIO

import docker
import pandas as pd
from skbio.stats.distance import mantel


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", action="store")
    parser.add_argument("-o", "--output_dir", action="store")
    results = parser.parse_args()
    return results.input_dir, results.output_dir


def main():
    log.basicConfig(level=log.INFO, filename="log_rarefaction.txt")
    log.info('Start')
    with open("results.txt", "w") as f:
        f.write("")

    input_dir, output_dir = get_arguments()
    shutil.rmtree(output_dir)

    read_numbers = [50, 100, 200, 500, 1000, 1500, 2000, 5000]

    for read_number in read_numbers:
        try:
            read_count = []
            os.mkdir(output_dir)
            log.info("reads_number {}".format(read_number))
            sampling_reads(input_dir, output_dir, read_number, sample_number=20)
            run_analysis("primary", "kmerprimary:1.0.0", output_dir)
            run_analysis("secondary", "kmersecondary:2.0.1", output_dir)
            mantel_array = get_coef_array(output_dir)
            log.info("Finish mantel analysis")

            for i in range(len(mantel_array)):
                read_count.append(read_number)

            with open("results.txt", "a") as f:
                f.write(str(mantel_array) + "\n")
                f.write(str(read_count) + "\n")
            shutil.rmtree(output_dir)
        except:
            log.info("Read number > reads in sample")

    log.info('Finished')


def sampling_reads(input_dir, output_dir, read_number, sample_number=10):
    log.info("Start sampling")
    files_array = os.listdir(input_dir)

    for index in range(sample_number):
        os.mkdir(os.path.join(output_dir, str(index)))

    for filename in files_array:
        with open(os.path.join(input_dir, filename), 'rb') as input_file:
            num_lines = sum([1 for line in input_file])
        total_records = int(num_lines / 4)

        output_files = []
        output_sequence_sets = []

        for i in range(sample_number):
            output_files.append(open(os.path.join(output_dir, str(i), filename[:-6]) + "_" + str(i) + ".fastq", "w"))
            output_sequence_sets.append(set(random.sample(range(total_records + 1), read_number)))

        record_number = 0
        with open(os.path.join(input_dir, filename), 'rb') as input_file:
            for line1 in input_file:
                line2 = next(input_file)
                line3 = next(input_file)
                line4 = next(input_file)
                for i, output in enumerate(output_files):
                    if record_number in output_sequence_sets[i]:
                        output.write(line1.decode())
                        output.write(line2.decode())
                        output.write(line3.decode())
                        output.write(line4.decode())
                record_number += 1
                if record_number % 100000 == 0:
                    print(str((record_number / total_records) * 100) + " % done")

        for output in output_files:
            output.close()
    log.info("Finish sampling")


def run_analysis(analysis, image, input_dir):
    log.info("Start {} analysis".format(analysis))
    for folder in os.listdir(input_dir):
        os.mkdir(os.path.join(input_dir, folder + analysis))
        abs_input_dir = os.path.abspath(os.path.join(input_dir, folder))
        abs_output_dir = os.path.abspath(os.path.join(input_dir, folder + analysis))
        docker_run(image, abs_input_dir, abs_output_dir)
        shutil.rmtree(abs_input_dir)
    log.info("Finish {} analysis".format(analysis))
    return True


def docker_run(image, input_dir, output_dir):
    client = docker.from_env()
    image = client.images.get(image)
    client.containers.run(image,
                          command='-i /docker_run/input -o /docker_run/output',
                          volumes={input_dir: {'bind': '/docker_run/input', 'mode': 'ro'},
                                   output_dir: {'bind': '/docker_run/output', 'mode': 'rw'}})
    return True


def get_coef_array(output_dir):
    dissimilarity_matrix_array = []
    coef_array = []

    log.info("Start mantel analysis")
    for folder in os.listdir(output_dir):
        with open(os.path.join(output_dir, folder, "basic_report_k_mer.json"), 'r') as json_file:
            kmer_json_str = json.load(json_file)
            kmer_json_dict = json.loads(kmer_json_str)
            kmer_link_content = kmer_json_dict["mds_kmer"]["kmer_dissim_matrix"]["link_content"]
            dissimilarity_matrix_array.append(get_data_frame_from_bytes(kmer_link_content))
    for i in range(len(os.listdir(output_dir))):
        for j in range(i + 1, len(os.listdir(output_dir))):
            coef, p_value, size = mantel(x=dissimilarity_matrix_array[i], y=dissimilarity_matrix_array[j],
                                         method='spearman', permutations=1000)
            coef_array.append(coef)
    return coef_array


def get_data_frame_from_bytes(string):
    df_string = base64.b64decode(string.encode("utf-8")).decode("utf-8")
    buffer = StringIO()
    buffer.write(df_string)
    buffer.seek(0)
    data_frame = pd.read_csv(buffer, index_col=0)
    return data_frame


if __name__ == "__main__":
    main()

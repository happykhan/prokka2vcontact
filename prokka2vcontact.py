#!/usr/bin/env python

# Copyright 2019 Nabil-Fareed Alikhan Licensed under the
#     Educational Community License, Version 2.0 (the "License"); you may
#     not use this file except in compliance with the License. You may
#     obtain a copy of the License at
#
#      http://www.osedu.org/licenses/ECL-2.0
#
#     Unless required by applicable law or agreed to in writing,
#     software distributed under the License is distributed on an "AS IS"
#     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
#     or implied. See the License for the specific language governing
#     permissions and limitations under the License.


"""
Converts Prokka output into a form for vContact2

### CHANGE LOG ###
2019-06-25 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build
2019-011-14 Evelien Adriaenssens
    * Clean-up
"""
import sys
import re
import os
import traceback
import argparse
import time
import logging

__licence__ = "ECL-2.0"
__author__ = "Nabil-Fareed Alikhan"
__author_email__ = "nabil@happykhan.com"
__version__ = "0.1.0"

epi = "Licence: " + __licence__ + " by " + __author__ + " <" + __author_email__ + ">"

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger()

def main():
    global args
    if args.verbose:
        logger.setLevel(logging.INFO)
    logger.info("Running...")
    faa_path = None
    table_path = None
    mapping_dict = {}
    for prok_file in os.listdir(args.prokka_dir):
        # TODO: Multiple files?
        if prok_file.endswith(".faa"):
            faa_path = os.path.join(args.prokka_dir, prok_file)
        if prok_file.endswith(".tbl"):
            table_path = os.path.join(args.prokka_dir, prok_file)
    if faa_path and table_path:
        logger.info("Extracting contig names from table file")
        count = 0
        with open(table_path) as table_f:
            contig_name = None
            locus_tag = None
            product_name = None
            for line in table_f.readlines():
                contig_match = re.search("^>Feature (.+)", line)
                if contig_match:
                    contig_name = contig_match.group(1).strip()
                product_match = re.search("\s+product\s+(.+)", line)
                if product_match:
                    product_name = product_match.group(1).strip()
                locus_match = re.search("\s+locus_tag\s+(.+)", line)
                if locus_match:
                    locus_tag = locus_match.group(1).strip()
                CDS_match = re.search("^\d+\s+\d+\s+\S+", line)
                if CDS_match:
                    if contig_name and locus_tag and product_name:
                        mapping_dict[locus_tag] = [contig_name, product_name]
                        real_CDS_match = re.search("^\d+\s+\d+\s+CDS", line)
                        if real_CDS_match:
                            count += 1
                    locus_tag = None
                    product_name = None
                    if count % 20000 == 0:
                        logger.info('Read %d CDS' % count)
            if contig_name and locus_tag and product_name:
                mapping_dict[locus_tag] = [contig_name, product_name]
        logger.info("Located %d CDS" % len(mapping_dict))
        logger.info("Reformatting faa")
        faa_out_path = args.output + ".faa"
        if args.csv:
            csv_out_path = args.output + ".csv"
        else:
            csv_out_path = args.output + ".tsv"

        with open(faa_path) as faa_f:
            faa_out = open(faa_out_path, "w")
            csv_out = open(csv_out_path, "w")
            if args.csv:
                csv_out.write("protein_id,contig_id,keywords\n")
            else:
                csv_out.write("protein_id\tcontig_id\tkeywords\n")
            faa_count = 0
            for line in faa_f.readlines():
                faa_header_match = re.search("^>(\S+)\s*.*", line)
                if faa_header_match:
                    faa_locus_tag = faa_header_match.group(1)
                    if mapping_dict.get(faa_locus_tag):
                        if faa_count % 20000 == 0:
                            logger.info('Written %d CDS' % faa_count)
                        faa_count += 1
                        faa_out.write(
                            ">%s %s [%s]\n"
                            % (
                                faa_locus_tag,
                                mapping_dict[faa_locus_tag][1],
                                mapping_dict[faa_locus_tag][0],
                            )
                        )
                        if args.csv:
                            csv_out.write(
                                "%s,%s,%s\n"
                                % (
                                    faa_locus_tag,
                                    mapping_dict[faa_locus_tag][0],
                                    mapping_dict[faa_locus_tag][1],
                                    )
                                )
                        else:
                            csv_out.write(
                                "%s\t%s\t%s\n"
                                % (
                                    faa_locus_tag,
                                    mapping_dict[faa_locus_tag][0],
                                    mapping_dict[faa_locus_tag][1],
                                    )
                                )
                        mapping_dict[faa_locus_tag].append('Done')
                    else:
                        logger.error('Locus tag %s is not found in prokka table file ' % faa_locus_tag)
                else:
                    faa_out.write(line)
        csv_out.close()
        faa_out.close()
        if faa_count != count:
            logger.error('Number of read CDS does not equal number of written, CDS include:')
            logger.error(','.join([x for x in mapping_dict if len(mapping_dict[x]) == 2]))
        logger.info("faa output is : %s" % os.path.abspath(faa_out_path))
        logger.info("csv output is : %s" % os.path.abspath(csv_out_path))
        logger.info('Done.')
    else:
        logger.error(
            "faa file or table file could not be found in %s . Exiting"
            % args.prokka_dir
        )


if __name__ == "__main__":
    try:
        start_time = time.time()
        desc = __doc__.split("\n\n")[0].strip()
        parser = argparse.ArgumentParser(description=desc, epilog=epi)
        parser.add_argument(
            "-v", "--verbose", action="store_true", default=False, help="verbose output"
        )
        parser.add_argument(
            "--version", action="version", version="%(prog)s " + __version__
        )
        parser.add_argument(
            "-o", "--output", action="store", help="output prefix", default="vc2"
        )
        parser.add_argument(
            "-c", "--csv", action="store_true", help="generate as csv ", default=False
        )
        parser.add_argument(
            "prokka_dir", action="store", help="Path to directory of Prokka output"
        )
        args = parser.parse_args()
        if args.verbose:
            print("Executing @ " + time.asctime())
        main()
        if args.verbose:
            print("Ended @ " + time.asctime())
        sys.exit(0)
    except KeyboardInterrupt:  # Ctrl-C
        raise
    except SystemExit:  # sys.exit()
        raise
    except Exception:
        print("ERROR, UNEXPECTED EXCEPTION")
        print(traceback.print_exc())
        os._exit(1)

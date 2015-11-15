#!/usr/bin/env python2.7

from bioblend.galaxy import GalaxyInstance
import requests
import datetime
import argparse
requests.packages.urllib3.disable_warnings()


def parse_args():
    args = argparse.ArgumentParser(description="Rename history items using a tabular file." +"\n" +
                                   "Example usage: python rename_hist_items.py -url misssissippi.snv.jussieu.fr \
                                   -key $your_api_key -hid $your_history_id -table $your_tabular_file \n \
                                   See test-data/sample_table.tab for an example file.")
    args.add_argument("-url", "--galaxy_url", required=True, help="url of galaxy instance")
    args.add_argument("-key", "--api_key", required=True, help="api key for galaxy instance" )
    args.add_argument("-hid", "--history_id", required=True, help="History id of hitory containing files to be renamed")
    args.add_argument("-table", "--rename_table", required=True, type=file, 
                      help="tab-seperated file with first column current filename,\
                      and second column the desired name")
    return args.parse_args()

    
def return_datetime(string_representation):
    """
    returns current time, to find last modified history.
    Currently ununsed, may be used in the future.
    """
    date, time = string_representation.split('T')
    return datetime.datetime.strptime(date + ' ' + time, "%Y-%m-%d %H:%M:%S.%f")


def get_rename_list(rename_table):
    return [(line.split('\t')[0],line.split('\t')[1].strip()) for line in rename_table]


def get_instance(url, api_key):
    return GalaxyInstance(url, api_key)


def get_name_id_d(gi, hid):
    return {dataset[u'name']:dataset[u'id'] for dataset in gi.histories.show_history(hid, contents=True)}


def update_names(gi, hid, rename_list, name_id_d ):
    for old_name, new_name in rename_list:
        dataset_id = name_id_d[old_name]
        gi.histories.update_dataset(history_id=hid, dataset_id=dataset_id, name=new_name)

        
def main():
    args = parse_args()
    hid = args.history_id
    rename_list = get_rename_list(args.rename_table)
    gi = get_instance(args.galaxy_url, args.api_key)
    name_id_d = get_name_id_d(gi, hid)
    rval = update_names(gi, hid, rename_list, name_id_d)

    
if __name__ == "__main__":
    main()

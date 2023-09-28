from collections import defaultdict
import os
import sys # for the command-line params
import time
import xml.etree.ElementTree as ET

import requests

import config
import db

connection = db.Connection()

def _record_data(data):
    """Parses a response from the efetch endpoint that has info about
    all the samples in the query."""
    tosave = None
    notfound = 0
    for package in data.findall('EXPERIMENT_PACKAGE'):
        sample = None

        for x in package.iter('SAMPLE'):
            if 'accession' in x.attrib.keys():
                sample = x.attrib['accession']
                break
        if sample is None:
            print("no sample found, moving on------------------------")
            continue

        for x in package.iter('INSTRUMENT_MODEL'):
            tosave = x.text

        if tosave is None:
            print("No instrument found, moving on--------------------")
            notfound += 1
            continue

        with connection.db.cursor() as cursor:
            cursor.execute("""
                UPDATE samples
                SET instrument=%s
                WHERE srs=%s
            """, (tosave, sample))
    return notfound

def find_instruments(count, per_query=50):
    """
    Queries the NCBI eUtils API to use sample IDs ("SRS" codes)
    to get information about runs ("SRR" codes) that can then
    be downloaded as FASTQ files.

    Inputs:
        - count: int. The upper limit for how many entries to search in total.
        - per_query: int. The number of entries to request in each web request
    """

    todo = connection.read("""
        SELECT srs FROM samples
        WHERE instrument IS NULL
        AND taxon='txid408170'
        LIMIT %s""", (count,))

    todo = [x[0] for x in todo] # each ID is nested inside a tuple of length 1
    print(f'Found {len(todo)} samples to process')
    cursor = 0
    notfound = 0
    while cursor < len(todo):
        if cursor > 0 and cursor % (per_query * 100) == 0:
            time.sleep(10) # pause for a bit every 100 requests
        if cursor % 1000 == 0:
            print(f'COMPLETE: {cursor}')

        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?tool={config.tool}&email={config.email}&db=sra&usehistory=y&term='
        for x in range(0, per_query):
            url += f'{todo[cursor]}[accn] or '
            cursor += 1
            if cursor == len(todo): break # in case the total isn't a multiple of "per_query"
        url = url[:-4] # trim off trailing " or "
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)
        try:
            r = requests.get(url)
        except:
            print('ERROR: Error sending request for webenv data. Skipping.')
            time.sleep(3)
            continue

        try:
            tree = ET.fromstring(r.text)
        except:
            print('ERROR: Couldnt parse response retrieving webenv data. Skipping.')
            time.sleep(3)
            continue

        webenv = tree.find('WebEnv')
        if webenv is None:
            print('\n---------\n')
            print(r.text)
            print("WARNING: Got response without a 'webenv' field. Moving on.")
            print('\n---\n')
            time.sleep(10)
            continue
        time.sleep(1)
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool={config.tool}&email={config.email}&db=sra&query_key=1&WebEnv={webenv.text}'
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)

        r = requests.get(url)
        try:
            tree = ET.fromstring(r.text)
        except ET.ParseError:
            print("WARNING: Misformed response from call to eFetch. Skipping.")
            time.sleep(10)
            continue
        notfound += _record_data(tree)
        #time.sleep(1)
    print(f'\n\n\nTOTAL MISSING INSTRUMENTS: {notfound} out of {count}')


if __name__ == "__main__":
    todo = 2000 if len(sys.argv) < 2 else sys.argv[1]

    find_instruments(todo, per_query=80)

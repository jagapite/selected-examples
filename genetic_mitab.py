#!/usr/bin/env python3
"""Script:: genetic_mitab.

Synopsis:
    Script that extracts FlyBase genetic interaction data from a reporting build 
    and writes in PSI MITAB format for export to Alliance.

Author:
    Julie Agapite jagapite@morgan.harvard.edu

"""

import psycopg2
import argparse
import csv
import re
import logging

# Set up the parser to handle cmd line arguments
parser = argparse.ArgumentParser(description='Query Chado.')
parser.add_argument('-s', '--pgserver', help='Postgres server', required=True)
parser.add_argument('-d', '--database', help='Reporting database', required=True)
parser.add_argument('-u', '--username', help='Postgres User name', required=True)
parser.add_argument('-p', '--password', help='Postgres password', required=True)
args = parser.parse_args()
server = args.pgserver
database = args.database
username = args.username
password = args.password

# Script must be run on a reporting db.
if not re.search('reporting$', database):
    raise argparse.ArgumentTypeError('Use of reporting database expected.')

# Generate filenames based on database used.
filename = "genetic_interactions_mitab_%s.tsv" % (database)
logfile = "genetic_interactions_mitab_%s.log" % (database)

# Set up logger.
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
file_handler = logging.FileHandler(logfile, mode='a')
formatter = logging.Formatter('%(asctime)s : %(levelname)s : Line No %(lineno)d : %(message)s')
file_handler.setFormatter(formatter)
log.addHandler(file_handler)

# Set up connection to db and query for genetic interactions.
conn_string = "host='%s' dbname='%s' user='%s' password='%s'" % (server, database, username, password)
conn = psycopg2.connect(conn_string)
main()
log.info('Genetic Interactions file is complete.')
conn.close()

def main():
    """Main function that gets all genetic interactions and pertinent info."""

    log.info('Starting query for interactions')

    # Query db and return a list of all interactions in form of tuples. [(fp_id, value, FBrf), (fp_id, value, FBrf), etc.]
    get_ints = ('SELECT DISTINCT fp.featureprop_id, fp.value, p.uniquename '
                'FROM featureprop fp, cvterm cvt, featureprop_pub fpp, pub p '
                'WHERE fp.type_id = cvt.cvterm_id '
                'AND cvt.name = \'derived_interaction_starting_gene\' '
                'AND fp.featureprop_id = fpp.featureprop_id AND fpp.pub_id = p.pub_id '
                'AND p.is_obsolete = \'f\' ORDER BY fp.value')
    int_tuples = connect(get_ints, 'no_query', conn)

    log.info("Featureprops found: %s", len(int_tuples))

    # Set up variables to hold interactions skipped and interactions for which alternate IDs could not be found.
    skipped = 0
    no_alt_id = set()
    no_pmid = set()
    no_related_FBrf = set()

    print_column_headers()

    for ixn in int_tuples:
        # Get info and print if interaction involves only Dmel genes (no \species) and is not a complex interaction with multiple genes affected (no @|@ or ::).
        if '@|@' in ixn[1] or '::' in ixn[1] or re.search(r'[A-Z][a-z]{3}\\', ixn[1]):
            skipped += 1
        else:
            # Make query tuple with FBrf to query db for pub info.
            if 'FBrf' in ixn[2]:
                pubid9 = 'flybase:' + ixn[2]
                pub_query = (ixn[2],)

                # Get author for column 8.
                get_author = ('SELECT DISTINCT pa.surname, pa.givennames, p.pyear '
                              'FROM pubauthor pa, pub p '
                              'WHERE p.pub_id = pa.pub_id AND pa.rank = 1 AND p.uniquename = %s')
                author = connect(get_author, pub_query, conn)
                if author:
                    if author[0][0] is not None:
                        last = author[0][0]
                    else:
                        last = ''
                    if author[0][1] is not None:
                        first = author[0][1]
                    else:
                        first = ''
                    if author[0][2] is not None:
                        year = author[0][2]
                    else:
                        year = ''
                    ref8 = last + " " + first + " (" + year + ")"
                else:
                    ref8 = '-'

                # Get pmid for column 9.
                # Use FBrf to query for pubmed ID.
                pmid = get_pubmed(pub_query, conn)
                if pmid:
                    pubid9 = pubid9 + '|pubmed:' + pmid
                else:
                    # Some FBrfs are a pub type that is not issued a pubmed ID, e.g. supplementary material. Look for related FBrfs that do have a pubmed ID.
                    get_related = ('SELECT DISTINCT p2.uniquename '
                                   'FROM pub p1, pub_relationship pr, pub p2 '
                                   'WHERE p1.pub_id = pr.object_id '
                                   'AND pr.subject_id = p2.pub_id '
                                   'AND p2.is_obsolete = \'f\' AND p1.uniquename = %s '
                                   'AND p2.uniquename LIKE \'FBrf%%\'')
                    related = connect(get_related, pub_query, conn)
                    if len(related) == 1:
                        pub_query = related

                    # Otherwise try another query to find FBrf of a related paper and use that to get pmid.
                    else:
                        get_also = ('SELECT DISTINCT p2.uniquename '
                                    'FROM pub p1, pub_relationship pr, pub p2 '
                                    'WHERE p1.pub_id = pr.subject_id '
                                    'AND pr.object_id = p2.pub_id '
                                    'AND p2.is_obsolete = \'f\' AND p1.uniquename = %s '
                                    'AND p2.uniquename LIKE \'FBrf%%\'')
                        also = connect(get_also, pub_query, conn)
                        if len(also) == 1:
                            pub_query = also
                        else:
                            no_related_FBrf.add(ixn[2])

                    pmid = get_pubmed(pub_query, conn)
                    if pmid:
                        pubid9 = pubid9 + '|pubmed:' + pmid
                    else:
                        no_pmid.add(ixn[2])
            else:
                # If interaction is not associated with an FBrf, it is unattributed and the pub and author fields are left blank.
                pubid9 = '-'
                ref8 = '-'

            # Construct dictionary to associate FBgns of generic genes to Entrez ID.
            # FB has generic genes that we use for curation of nearly identical gene
            # families when the specific family member cannot be determined. In order
            # to include these interactions, one member of the gene family was chosen
            # as representative.
            gene_map = {}
            gene_map['FBgn0000002'] = '3771903'
            gene_map['FBgn0001195'] = '3772715'
            gene_map['FBgn0001196'] = '318855'
            gene_map['FBgn0001198'] = '326273'
            gene_map['FBgn0001199'] = '318847'
            gene_map['FBgn0001200'] = '318846'
            gene_map['FBgn0003523'] = '2768872'
            gene_map['FBgn0061471'] = '5740577'
            gene_map['FBgn0061474'] = '26067172'
            gene_map['FBgn0061475'] = '5740812'
            gene_map['FBgn0065042'] = '26067164'

            # Extract interaction participants from featureprop statement.
            parts = re.search('^@(FBgn[0-9]{7}):(.+)@ (.+) @(FBgn[0-9]{7}):(.+)@', ixn[1])

            # Primary ID for column 1 is FBgn.
            A = parts.group(1)
            A_id1 = 'flybase:' + A

            # If odd characters in gene name for column 5, place name in quotes.
            if '(' in parts.group(2) or ':' in parts.group(2):
                geneA = '"' + parts.group(2) + '"'
            else:
                geneA = parts.group(2)
            A_name5 = 'flybase:' + geneA + '(gene name)'

            # Make tuple to use in queries based on gene ID (FBgn).
            gid = (A,)

            # Get alt IDs for column 3.
            CG = get_CG_id(gid, conn)
            if A in gene_map:
                Entrez = gene_map[A]
            else:
                Entrez = get_Entrez_id(gid, conn)
            if CG and Entrez:
                A_altids3 = 'flybase:' + CG[0][0] + '|entrez gene/locuslink:' + Entrez
            elif Entrez and not CG:
                A_altids3 = 'entrez gene/locuslink:' + Entrez
            elif CG and not Entrez:
                A_altids3 = 'flybase:' + CG[0][0]
            else:
                A_altids3 = '-'
                no_alt_id.add(A)

            # Get taxid for column 10.
            tax = get_taxid(gid, conn)
            if tax is not None:
                A_tax10 = 'taxid:' + tax[0] + '(' + tax[1] + ' ' + tax[2] + ')'
            else:
                A_tax10 = '-'

            # Get same information for interactor B.
            B = parts.group(4)
            B_id2 = 'flybase:' + B

            if '(' in parts.group(5) or ':' in parts.group(5):
                geneB = '"' + parts.group(5) + '"'
            else:
                geneB = parts.group(5)
            B_name6 = 'flybase:' + geneB + '(gene name)'

            gid = (B,)

            CG = get_CG_id(gid, conn)
            if B in gene_map:
                Entrez = gene_map[B]
            else:
                Entrez = get_Entrez_id(gid, conn)
            if CG and Entrez:
                B_altids4 = 'flybase:' + CG[0][0] + '|entrez gene/locuslink:' + Entrez
            elif Entrez and not CG:
                B_altids4 = 'entrez gene/locuslink:' + Entrez
            elif CG and not Entrez:
                B_altids4 = 'flybase:' + CG[0][0]
            else:
                B_altids4 = '-'
                no_alt_id.add(B)

            tax = get_taxid(gid, conn)
            if tax is not None:
                B_tax11 = 'taxid:' + tax[0] + '(' + tax[1] + ' ' + tax[2] + ')'
            else:
                B_tax11 = '-'

            # Get interaction type and interactor roles. Role terms are pending incorporation into the psi-mi ontology.
            if 'suppressible' in parts.group(3):
                intype12 = 'psi-mi:"MI:0796"(suppression)'
                A_role19 = 'psi-mi:"MI:0582"(suppressed gene)'
                B_role20 = 'psi-mi:"MI:0581"(suppressor gene)'
            elif 'enhanceable' in parts.group(3):
                intype12 = 'psi-mi:"MI:1271"(enhancement)'
                A_role19 = 'psi-mi:"MI:2352"(enhanced gene)'
                B_role20 = 'psi-mi:"MI:2351"(enhancer gene)'
            else:
                intype12 = '-'
                A_role19 = '-'
                B_role20 = '-'

            # Print output row for this interaction.
            with open(filename, 'a') as csvfile:
                csv_writer = csv.writer(csvfile, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
                csv_writer.writerow([A_id1, B_id2, A_altids3, B_altids4,
                                     A_name5, B_name6, 'psi-mi:"MI:0254"(genetic interference)', ref8,
                                     pubid9, A_tax10, B_tax11, intype12, 'psi-mi:"MI:0478"(flybase)', '-', '-',
                                     '-', '-', '-', A_role19, B_role20, 'psi-mi:"MI:0250"(gene)',
                                     'psi-mi:"MI:0250"(gene)', '-', '-', '-', '-', '-', '-', '-',
                                     '-', '-', '-', '-', '-', '-', 'false', '-', '-', '-',
                                     '-', '-', '-'])

    log.info("Number of interactions skipped: %s", skipped)
    log.info("Number of genes with no alt ids: %s", len(no_alt_id))
    log.info("Total number of FBrfs with no pmid: %s", len(no_pmid))
    log.info("Number of FBrfs with no related FBrf: %s", len(no_related_FBrf))
    log.info("FBrfs with related FBrf but no pmid: %s", no_pmid.difference(no_related_FBrf))


def connect(sql, query, conn):
    """Query database and return results of query."""
    cursor = conn.cursor()
    if query == 'no_query':
        cursor.execute(sql)
    else:
        cursor.execute(sql, query)
    records = cursor.fetchall()
    cursor.close()
    return records


def print_column_headers():
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        csv_writer.writerow(['#ID(s) Interactor A', 'ID(s) Interactor B', 'Alt ID(s) Interactor A',
                             'Alt ID(s) Interactor B', 'Alias(es) Interactor A', 'Alias(es) Interactor B',
                             'Interaction Detection Method(s)', 'Publication 1st Author(s)', 'Publication ID(s)',
                             'Taxid Interactor A', 'Taxid Interactor B', 'Interaction Type(s)', 'Source Database(s)',
                             'Interaction Identifier(s)', 'Confidence Value(s)', 'Expansion Method(s)',
                             'Biological Role(s) Interactor A', 'Biological Role(s) Interactor B',
                             'Experimental Role(s) Interactor A', 'Experimental Role(s) Interactor B',
                             'Type(s) Interactor A', 'Type(s) Interactor B', 'Xref(s) Interactor A',
                             'Xref(s) Interactor B', 'Interaction Xref(s)', 'Annotation(s) Interactor A',
                             'Annotation(s) Interactor B', 'Interaction Annotation(s)', 'Host Organism(s)',
                             'Interaction Parameters', 'Creation Date', 'Update Date', 'Checksum Interactor A',
                             'Checksum Interactor B', 'Interaction Checksum', 'Negative', 'Feature(s) Interactor A',
                             'Feature(s) Interactor B', 'Stoichiometry Interactor A', 'Stoichiometry Interactor B',
                             'Identification Method(s) Participant A', 'Identification Method(s) Participant B'])


def get_pubmed(pub_query, conn):
    """Return pubmed id for the input FBrf."""
    get_pm = ('SELECT DISTINCT dx.accession '
              'FROM pub p, pub_dbxref pd, dbxref dx, db '
              'WHERE p.pub_id = pd.pub_id '
              'AND pd.dbxref_id = dx.dbxref_id AND dx.db_id = db.db_id '
              'AND db.name = \'pubmed\' AND pd.is_current = TRUE '
              'AND p.uniquename =  %s')
    pm = connect(get_pm, pub_query, conn)
    if pm:
        return(pm[0][0])
    else:
        return(None)


def get_CG_id(gid, conn):
    """Return CG ID for input FBgn."""
    get_CG = ('SELECT DISTINCT dx.accession '
              'FROM feature f, feature_dbxref fd, db, dbxref dx '
              'WHERE f.feature_id = fd.feature_id AND fd.dbxref_id = dx.dbxref_id '
              'AND dx.db_id = db.db_id AND db.name = \'FlyBase Annotation IDs\' AND '
              'dx.accession NOT LIKE \'%%-%%\' AND fd.is_current = \'t\' '
              'AND f.uniquename = %s')
    CG_id = connect(get_CG, gid, conn)
    return(CG_id)


def get_Entrez_id(gid, conn):
    """Return Entrez ID for input FBgn."""
    get_Entrez = ('SELECT DISTINCT dx.accession '
                  'FROM feature f, feature_dbxref fd, db, dbxref dx '
                  'WHERE f.feature_id = fd.feature_id AND fd.dbxref_id = dx.dbxref_id '
                  'AND dx.db_id = db.db_id AND db.name = \'EntrezGene\' AND '
                  'fd.is_current = \'t\' AND f.uniquename = %s')
    Entrez_id = connect(get_Entrez, gid, conn)
    if Entrez_id:
        id = Entrez_id[0][0]
    else:
        id = None
    return(id)


def get_taxid(gid, conn):
    """Return NCBI taxid, genus, species for input feature ID."""
    taxid = ('SELECT DISTINCT dx.accession, o.genus, o.species '
             'FROM feature g, organism o, organism_dbxref od, dbxref dx, db '
             'WHERE g.organism_id = o.organism_id AND o.organism_id = od.organism_id '
             'AND od.dbxref_id =dx.dbxref_id AND dx.db_id = db.db_id '
             'AND db.name = \'NCBITaxon\' '
             'AND g.uniquename = %s')
    tid = connect(taxid, gid, conn)
    if len(tid) > 0:
        return(tid[0][0], tid[0][1], tid[0][2])
    else:
        return(None)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import psycopg2
import argparse
import csv
import re


def connect(sql, query, conn):
    cursor = conn.cursor() # Return the cursor and use it to perform queries.

    # Execute the query.
    if query == 'no_query': # If we're just running an SQL query without a variable, only execute the sql
        cursor.execute(sql)
    else:
        cursor.execute(sql, query) # If we have a variable (e.g. an FBgn) be sure to include it in the execute command.
    
    records = cursor.fetchall() # Grab the results.
    cursor.close() # Close the cursor.
    return records # Return a list of tuples.


def get_pubmed(pub_query,conn):
    # Returns pubmed id for the input FBrf

    get_pm = ('SELECT DISTINCT dx.accession '
                        'FROM pub p, pub_dbxref pd, dbxref dx, db '
                        'WHERE p.pub_id = pd.pub_id '
                        'AND pd.dbxref_id = dx.dbxref_id AND dx.db_id = db.db_id '
                        'AND db.name = \'pubmed\' AND pd.is_current = TRUE AND p.uniquename =  %s')
      
    pm = connect(get_pm, pub_query, conn)
    if pm:
        return(pm[0][0])
    else:
        return(None)


def get_CG_id(gid, conn):
    # Returns CG ID for input FBgn

    get_CG = ('SELECT DISTINCT dx.accession '
                'FROM feature f, feature_dbxref fd, db, dbxref dx '
                'WHERE f.feature_id = fd.feature_id AND fd.dbxref_id = dx.dbxref_id '
                'AND dx.db_id = db.db_id AND db.name = \'FlyBase Annotation IDs\' AND '
                'dx.accession NOT LIKE \'%%-%%\' AND fd.is_current = \'t\' AND f.uniquename = %s')
    CG_id = connect(get_CG,gid,conn)
    return(CG_id)


def get_Entrez_id(gid,conn):
    # Returns Entrez ID for input FBgn

    get_Entrez = ('SELECT DISTINCT dx.accession '
                    'FROM feature f, feature_dbxref fd, db, dbxref dx '
                    'WHERE f.feature_id = fd.feature_id AND fd.dbxref_id = dx.dbxref_id '
                    'AND dx.db_id = db.db_id AND db.name = \'EntrezGene\' AND '
                    'fd.is_current = \'t\' AND f.uniquename = %s')
    Entrez_id = connect(get_Entrez,gid,conn)
    if Entrez_id:
        id = Entrez_id[0][0]
    else:
        id = None
    return(id)


def get_taxid(xid,conn):
    # Returns NCBI taxid, genus, species for input feature ID

    taxid = ('SELECT DISTINCT dx.accession, o.genus, o.species '
            'FROM feature x, organism o, organism_dbxref od, dbxref dx, db '
            'WHERE x.organism_id = o.organism_id AND o.organism_id = od.organism_id '
            'AND od.dbxref_id =dx.dbxref_id AND dx.db_id = db.db_id AND db.name = \'NCBITaxon\' '
            'AND x.uniquename = %s')
    tid = connect(taxid,xid,conn)
    if len(tid) > 0:
        return(tid[0][0],tid[0][1],tid[0][2])
    else:
        return(None)


def get_role(xint,conn):
    # Returns participant role for input interaction and participant

    role = ('SELECT DISTINCT dx.accession, cvt.name, cvt.cvterm_id '
                    'FROM interaction i, feature_interaction fi, feature f, cvterm cvt, dbxref dx '
                    'WHERE f.feature_id = fi.feature_id AND fi.interaction_id = i.interaction_id '
                    'AND fi.role_id = cvt.cvterm_id AND cvt.dbxref_id = dx.dbxref_id '
                    'AND f.uniquename = %s AND i.uniquename = %s')
    part_role = connect(role,xint,conn)
    return(part_role)


def get_subregions(xint,conn): 
    # Returns subregions associated with input interaction and participant

    subregions = ('SELECT DISTINCT cvt.name, fip.value, f.name '
                    'FROM interaction i, feature_interaction fi, feature_interactionprop fip, ' 
                    'feature f, cvterm cvt, cvterm cvt2, feature_relationship fr, feature f2 '
                    'WHERE f.feature_id = fi.feature_id AND fi.interaction_id = i.interaction_id '
                    'AND fi.feature_interaction_id = fip.feature_interaction_id '
                    'AND fi.role_id = cvt.cvterm_id '
                    'AND fip.type_id = cvt2.cvterm_id AND '
                    'cvt2.name = \'subpart_info\' AND f.feature_id = fr.subject_id '
                    'AND f2.feature_id = fr.object_id AND f.is_obsolete = \'f\' AND '
                    'f2.uniquename = %s AND i.uniquename = %s')
    subs = connect(subregions,xint,conn)
    return(subs)


def get_isoforms(xint,conn):
    # Returns isoform associated with an interaction and participant

    isoforms = ('SELECT DISTINCT f.name '
                    'FROM interaction i, feature_interaction fi, feature_interactionprop fip, ' 
                    'feature f, cvterm cvt, cvterm cvt2, feature_relationship fr, feature f2 '
                    'WHERE f.feature_id = fi.feature_id AND fi.interaction_id = i.interaction_id '
                    'AND fi.feature_interaction_id = fip.feature_interaction_id '
                    'AND fi.role_id = cvt.cvterm_id '
                    'AND fip.type_id = cvt2.cvterm_id AND '
                    'cvt2.name = \'interacting isoform\' AND f.feature_id = fr.subject_id '
                    'AND f2.feature_id = fr.object_id AND f.is_obsolete = \'f\' AND '
                    'f2.uniquename = %s AND i.uniquename = %s')
    isos = connect(isoforms,xint,conn)
    return(isos)


def get_tag_info(xint,conn): 
    # Returns tag/experimental feature/notes field associated with an interaction and participant

    get_tags = ('SELECT DISTINCT fip2.value '
            'FROM interaction i, feature_interaction fi, feature_interactionprop fip, '
            'feature f, cvterm cvt, feature_interactionprop fip2, cvterm cvt2 '
            'WHERE f.feature_id = fi.feature_id AND fi.interaction_id = i.interaction_id '
            'AND fi.feature_interaction_id = fip.feature_interaction_id '
            'AND fip.type_id = cvt.cvterm_id AND cvt.name = \'participating feature\' '
            'AND fi.feature_interaction_id = fip2.feature_interaction_id AND fip2.type_id = cvt2.cvterm_id '
            'AND cvt2.name = \'comment\' AND f.uniquename = %s AND i.uniquename = %s')
    tags = connect(get_tags,xint,conn)
    return(tags)


def get_child_ids(id,conn):
    # Returns list of all child terms of input cvterm id

    child_ids = ('WITH RECURSIVE children AS '
                    '(SELECT subject_id '
                    'FROM cvterm_relationship '
                    'WHERE object_id = %s '
                    'UNION '
                    'SELECT cr.subject_id '
                    'FROM cvterm_relationship cr '
                    'INNER JOIN children ch ON ch.subject_id = cr.object_id) '
                    'SELECT * FROM children')
    ids = connect(child_ids,id,conn)
    list_of_ids = []
    for item in ids:
        list_of_ids.append(item[0])
    return(list_of_ids)


def  get_comments(qint,conn):
    # Get all other comments except source and internal notes

    comms = ('SELECT DISTINCT ip.value '
                    'FROM interaction i, interactionprop ip, cvterm cvt '
                    'WHERE i.interaction_id = ip.interaction_id AND ip.type_id = cvt.cvterm_id '
                    'AND cvt.is_obsolete=0 AND cvt.name != \'comments on source\' '
                    'AND cvt.name != \'internalnotes\' AND i.uniquename = %s')
    comnts = connect(comms, qint, conn)
    return(comnts)


def whittle(xlist,extra_terms):
    # Removes extra terms from input list
    for term in xlist:
        if term[1] in extra_terms:
            xlist.remove(term)   
    return(xlist)


def query_for_ints(filename, conn): 
    # Main function- gets all interactions and info

    # Query db and return a list of all interactions in form of tuples. [(int, FBig), (int, FBig), etc.]
    get_ints = ('SELECT DISTINCT i.uniquename, ig.uniquename '
                        'FROM interaction i, feature_interaction fi, interaction_group_feature_interaction igfi, interaction_group ig  '
                        'WHERE i.is_obsolete=\'f\' AND ig.is_obsolete = \'f\' AND i.interaction_id = fi.interaction_id '
                        'AND fi.feature_interaction_id = igfi.feature_interaction_id AND igfi.interaction_group_id = ig.interaction_group_id '
                        'ORDER BY ig.uniquename')
    int_tuples = connect(get_ints, 'no_query', conn)

    # Print column headers
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter = '\t')
        csv_writer.writerow(['#ID(s) Interactor A', 'ID(s) Interactor B','Alt ID(s) Interactor A', 'Alt ID(s) Interactor B',
            'Alias(es) Interactor A', 'Alias(es) Interactor B', 'Interaction Detection Method(s)',
            'Publication 1st Author(s)', 'Publication ID(s)', 'Taxid Interactor A', 'Taxid Interactor B',
            'Interaction Type(s)', 'Source Database(s)', 'Interaction Identifier(s)', 'Confidence Value(s)',
            'Expansion Method(s)', 'Biological Role(s) Interactor A', 'Biological Role(s) Interactor B',
            'Experimental Role(s) Interactor A', 'Experimental Role(s) Interactor B', 'Type(s) Interactor A',
            'Type(s) Interactor B', 'Xref(s) Interactor A', 'Xref(s) Interactor B', 'Interaction Xref(s)',
            'Annotation(s) Interactor A', 'Annotation(s) Interactor B', 'Interaction Annotation(s)', 'Host Organism(s)',
            'Interaction Parameters', 'Creation Date', 'Update Date', 'Checksum Interactor A', 'Checksum Interactor B',
            'Interaction Checksum', 'Negative', 'Feature(s) Interactor A', 'Feature(s) Interactor B', 'Stoichiometry Interactor A',
            'Stoichiometry Interactor B', 'Identification Method(s) Participant A', 'Identification Method(s) Participant B'])
    
    # Get a list of detection method ids that are children of 'interaction detection method'
    # to distinguish from participant detection methods since we curate both types into same field
    method_tuple = (103168,)
    int_method_list = get_child_ids(method_tuple,conn)

    # Get a list of experimental role ids that are children of 'experimental role'
    # to distinguish from biological roles since we curate both types into same field
    role_tuple = (103631,)
    role_list = get_child_ids(role_tuple,conn)

    # For each interaction, get data
    for int in int_tuples:

        # Place interaction ID in variable for column 14
        int_id14 = 'flybase:' + int[0]

        # Place FBig ID in variable for column 25
        int_xref25 = 'flybase:' + int[1]

        # Make a query tuple with interaction ID to get more data
        qint = (int[0],)

        # Get collection ID associated with interaction if there is one and add to column 25
        get_FBlc = ('SELECT DISTINCT l.uniquename '
                        'FROM library l, library_interaction li, interaction i '
                        'WHERE i.interaction_id = li.interaction_id AND li.library_id = l.library_id '
                        'AND i.is_obsolete = \'f\' AND l.is_obsolete = \'f\' AND i.uniquename = %s')
        FBlc = connect(get_FBlc, qint, conn)
        if FBlc:
            int_xref25 += '|flybase:' + FBlc[0][0]

        # Get author for column 8
        get_author = ('SELECT DISTINCT pa.surname, pa.givennames, p.pyear '
                        'FROM pubauthor pa, pub p, interaction_pub ip, interaction i '
                        'WHERE i.interaction_id = ip.interaction_id AND ip.pub_id= p.pub_id '
                        'AND p.pub_id = pa.pub_id AND pa.rank = 1 AND i.uniquename = %s')
        author = connect(get_author, qint, conn)
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

        # Get pub IDs, FBrf and pmid, for column 9
        # Get FBrf and pub type -all interactions should be associated with an FBrf
        get_FBrf = ('SELECT DISTINCT p.uniquename, cvt.name '
                        'FROM interaction i, interaction_pub ip, pub p, cvterm cvt '
                        'WHERE i.interaction_id = ip.interaction_id '
                        'AND ip.pub_id = p.pub_id AND p.type_id = cvt.cvterm_id AND '
                        'p.is_obsolete = \'f\' AND i.uniquename =  %s')
        FBrf = connect(get_FBrf, qint, conn)
        pubid9 = 'flybase:' + FBrf[0][0]

        # Use FBrf to query for pubmed ID
        # If associated FBrf is a personal communication, it will not have a pmid
        if FBrf[0][1] == 'personal communication to FlyBase':
            pub_query = ()

        # If associated FBrf is a paper and not e.g. suppl material, use FBrf to get pmid
        elif FBrf[0][1] == 'paper':
            pub_query = (FBrf[0][0],)

        # If associated FBrf is not itself a paper, get FBrf of related paper
        else:
            get_related = ('SELECT DISTINCT p2.uniquename '
                            'FROM pub p1, pub_relationship pr, pub p2, cvterm cvt '
                            'WHERE p1.pub_id = pr.object_id AND pr.subject_id = p2.pub_id '
                            'AND p2.type_id = cvt.cvterm_id AND p2.is_obsolete = \'f\' AND p1.uniquename = %s')
            rel_query = (FBrf[0][0],)
            related = connect(get_related,rel_query,conn)
            if related:
                pub_query = related

        # If FBrf for a paper is found, use it to get pmid
        if pub_query:
            pmid = get_pubmed(pub_query,conn)
            if pmid:
                pubid9 = pubid9 + '|pubmed:' + pmid
            
            # Otherwise try another query to find FBrf of a related paper and use that to get pmid
            else:
                get_also = ('SELECT DISTINCT p2.uniquename '
                            'FROM pub p1, pub_relationship pr, pub p2, cvterm cvt '
                            'WHERE p1.pub_id = pr.subject_id AND pr.object_id = p2.pub_id '
                            'AND pr.type_id = cvt.cvterm_id AND cvt.name = \'also_in\' '
                            'AND p2.is_obsolete = \'f\' AND p1.uniquename = %s')
                also = connect(get_also,pub_query,conn)
                if also:
                    pub_query = also
                    pmid = get_pubmed(pub_query,conn)
                    if pmid:
                        pubid9 = pubid9 + '|pubmed:' + pmid

        # Get interaction type for column 12
        get_int_type = ('SELECT dx.accession, cvt.name '
                        'FROM interaction i, cvterm cvt, dbxref dx, db '
                        'WHERE i.type_id = cvt.cvterm_id AND cvt.dbxref_id = dx.dbxref_id '
                        'AND dx.db_id = db.db_id AND db.name = \'MI\' AND i.uniquename = %s')
        int_type = connect(get_int_type, qint, conn)
        intype12 = 'psi-mi:"MI:' + int_type[0][0] + '"(' + int_type[0][1] + ')'

        # Get interaction annotations for column 28
        # First get comments on source
        get_source = ('SELECT DISTINCT ip.value '
                        'FROM interaction i, interactionprop ip, cvterm cvt '
                        'WHERE i.interaction_id = ip.interaction_id AND ip.type_id = cvt.cvterm_id '
                        'AND cvt.is_obsolete=0 AND cvt.name = \'comments on source\' '
                        'AND i.uniquename = %s')
        annots = connect(get_source, qint, conn)

        # If there are source annotations, start the string with the first
        if annots:
            annots28 = 'comment:"' + annots[0][0].replace('"','') + '"'

            # If more than one source annotation, add others to string
            if len(annots) > 1:
                for ann in annots[1:]:
                    annots28 += '|comment:"' + ann[0].replace('"','') + '"'

            # Then get comments and add to string
            comments = get_comments(qint,conn)
            if comments:
                for comm in comments:
                    if comm[0] is not None:
                        annots28 += '|comment:"' + comm[0].replace('"','') + '"'

        # If no source annotations, get comments and start string with comment
        else:
            comments = get_comments(qint,conn)
            if comments:
                annots28 = 'comment:"' + comments[0][0].replace('"','') + '"'
                if len(comments) > 1:
                    for comm in comments[1:]:
                        annots28 += '|comment:"' + comm[0].replace('"','') + '"'
            else:
                annots28 = '-'

        # Get interaction detection method for column 7
        get_methods = ('SELECT DISTINCT dx.accession, cvt.name, cvt.cvterm_id '
                        'FROM dbxref dx, cvterm cvt, interaction_cvterm ic, interaction i '
                        'WHERE dx.dbxref_id = cvt.dbxref_id AND cvt.cvterm_id = ic.cvterm_id AND '
                        'ic.interaction_id = i.interaction_id AND cvt.is_obsolete=0 AND i.uniquename = %s')
        methods = connect(get_methods,qint,conn)
        assays7 = ''

        # Start a list of methods attached to interaction that are in approprate cv branch
        ok_list = []
        for meth in methods:
            if meth[2] in int_method_list:
                ok_list.append(meth)

        # Make lists of terms to remove if there is more than one appropriate method term        
        term_list_1 = ['inferred by author', 'inferred by curator', 'competition binding', 'genetic interference', 'light microscopy',
                    'ultraviolet-visible spectroscopy', 'interologs mapping', 'luminiscence technology', 'experimental knowledge based']
        term_list_2 = ['fluorescence microscopy', 'nucleic acid uv cross-linking assay', 'cross-linking study', 'fluorescence technology',
                    'protein cross-linking with a bifunctional reagent', 'competition binding', 'confocal microscopy']

        # If there is only one method on the ok_list make a string
        if len(ok_list) == 1:
            assays7 = 'psi-mi:"MI:' + ok_list[0][0] + '"(' + ok_list[0][1] + ')'
        
        # If there is more than one method, remove terms from list one, then list two if still more than one
        elif len(ok_list) > 1:
            better_list = whittle(ok_list, term_list_1)
            if len(better_list) == 1:
                assays7 = 'psi-mi:"MI:' + better_list[0][0] + '"(' + better_list[0][1] + ')'
            elif len(better_list) > 1:
                best_list = whittle(better_list, term_list_2)
                if len(best_list) == 1:
                    assays7 = 'psi-mi:"MI:' + best_list[0][0] + '"(' + best_list[0][1] + ')'

                # If after removing extraneous terms there is still more than one method, choose one based on
                # hierarchy affinity > enzymatic > cosedimentation > scattering
                elif len(best_list) > 1:
                    for item in best_list:
                        if ('affinity' in item[1] or 'coimmunoprecipitation' in item[1] or 'pull down' in item[1]):
                            assays7 = 'psi-mi:"MI:' + item[0] + '"(' + item[1] + ')'
                            break                
                    if not assays7:
                        for item in best_list:
                            if 'enzymatic' in item[1] and not assays7:
                                assays7 = 'psi-mi:"MI:' + item[0] + '"(' + item[1] + ')'
                                break
                        if not assays7:
                            for item in best_list:
                                if 'cosedimentation' in item[1] and not assays7:
                                    assays7 = 'psi-mi:"MI:' + item[0] + '"(' + item[1] + ')'
                                    break
                            if not assays7:
                                for item in best_list:
                                    if 'scattering' in item[1] and not assays7:
                                        assays7 = 'psi-mi:"MI:' + item[0] + '"(' + item[1] + ')'
                                        break

                    # If a single method term has not yet been identified, choose based on interaction id
                    if not assays7:
                        if 'CH' in int[0]:
                            assays7 = 'psi-mi:"MI:0091"(chromatography technology)'
                        elif 'FT' in int[0]:
                            assays7 = 'psi-mi:"MI:0417"(footprinting)'
                        
                        # Otherwise use multiple terms
                        else:
                            assays7 = 'psi-mi:"MI:' + best_list[0][0] + '"(' + best_list[0][1] + ')'
                else:
                    assays7 = '-'
            else:
                assays7 = '-'
        else:
            assays7 = '-'

        # Get interaction participants
        get_parts = ('SELECT DISTINCT g.uniquename, g.name, x.uniquename, x.name '
                    'FROM feature g, feature_relationship fr, feature x, feature_interaction fi, '
                    'feature_interactionprop fip, interaction i, cvterm cvt '
                    'WHERE i.interaction_id = fi.interaction_id AND fi.feature_interaction_id = fip.feature_interaction_id '
                    'AND fip.type_id = cvt.cvterm_id AND fi.feature_id = x.feature_id '
                    'AND x.feature_id = fr.subject_id AND fr.object_id = g.feature_id AND fr.type_id = 59983 '
                    'AND cvt.name = \'participating feature\' AND g.is_obsolete = \'f\' AND x.is_obsolete = \'f\' '
                    'AND g.uniquename LIKE \'FBgn%%\' AND i.uniquename = %s ORDER BY g.uniquename')
        parts = connect(get_parts,qint,conn)

        # Construct dictionary to associate FBgns of generic genes to Entrez ID
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

        # Get information for interactor A
        # Primary ID for column 1 id FBgn
        A_id1 = 'flybase:' + parts[0][0]

        # If odd characters in gene name, place name in quotes
        if '(' in parts[0][1] or ':' in parts[0][1]:
            geneA = '"' + parts[0][1] + '"'
        else:
            geneA = parts[0][1]
        A_name5 = 'flybase:' + geneA + '(gene name)'

        # Assign molecule type based on feature name
        if '-XP' in parts[0][3]:
            A_type21 = 'psi-mi:"MI:0326"(protein)'
        elif '-XR' in parts[0][3]:
                A_type21 = 'psi-mi:"MI:0320"(ribonucleic acid)'
        else:
            A_type21 = 'psi-mi:"MI:0329"(unknown participant)'

        # Make tuple to use in queries based on gene ID (FBgn)
        gid = (parts[0][0],)

        # Get alt IDs for column 3
        CG = get_CG_id(gid,conn)
        if parts[0][0] in gene_map:
            Entrez = gene_map[parts[0][0]]
        else:
            Entrez = get_Entrez_id(gid,conn)

        if CG and Entrez:
            A_altids3 = 'flybase:' + CG[0][0] + '|entrez gene/locuslink:' + Entrez
        elif Entrez and not CG:
            A_altids3 = 'entrez gene/locuslink:' + Entrez
        elif CG and not Entrez:
            A_altids3 = 'flybase:' + CG[0][0]
        else:
            A_altids3 = '-'

        # Make tuple to use in queries based on feature ID
        xid = (parts[0][2],)

        # Get taxid for column 10
        tax = get_taxid(xid,conn)
        if tax is not None:
            A_tax10 = 'taxid:' + tax[0] + '("' + tax[1] + ' ' + tax[2] + '")' 
        else:
             A_tax10 = '-'

        # Make tuple for queries based on interaction ID and feature ID
        xint = (parts[0][2],int[0],)

        rol = get_role(xint,conn)
        if rol[0][2] in role_list:
            A_role19 =  'psi-mi:"MI:' + rol[0][0] + '"(' + rol[0][1] + ')'
            A_role17 = '-'
        else:
            A_role17 = 'psi-mi:"MI:' + rol[0][0] + '"(' + rol[0][1] + ')'
            A_role19 = 'psi-mi:"MI:0499"(unspecified role)'

    #    Comment out section to get subregions until data can be sufficiently standardized for mitab parser
    #    Get subregion/experimental feature info for A
    #    Asubs = get_subregions(xint,conn)
    #    if Asubs:
    #        if 'aa ' in Asubs[0][1] or 'nt ' in Asubs[0][1]:
    #            Asubregion = Asubs[0][1][3:]
    #            if 'aa ' in Asubregion:
    #                Asubregion = Asubregion.replace(" and aa ", ",")
    #        else:
    #            Asubregion = Asubs[0][1]
    #        regexp = re.compile(r'[0-9]+,[0-9]+,')
    #        if ('-' in Asubregion and regexp.search(Asubregion)) or ('-' not in Asubregion and ',' in Asubregion):
    #            components = Asubregion.split(",")
    #            Asubregion = components[0] + "-" + components[len(components)-1]
    #        elif '-' not in Asubregion and ',' not in Asubregion:
    #            Asubregion = Asubregion + "-" + Asubregion
    #        if "(" in Asubs[0][2] or "-" in Asubs[0][2]:
    #            Atext = '"' + Asubs[0][2] + '"'
    #        else:
    #            Atext = Asubs[0][2]
    #        A_feature37 = Asubs[0][0] + ':' + Asubregion + '(' + Atext + ')'
    #        if len(Asubs) > 1:
    #            for item in Asubs[1:]:
    #                if 'aa ' in item[1] or 'nt ' in item[1]:
    #                    Asub = item[1][3:]
    #                    if 'aa ' in Asub:
    #                        Asub = Asub.replace(" and aa ", ",")
    #                else:
    #                    Asub = item[1]
    #                if ('-' in Asub and regexp.search(Asub)) or ('-' not in Asub and ',' in Asub):
    #                    components = Asub.split(",")
    #                    Asub = components[0] + "-" + components[len(components)-1]
    #                elif '-' not in Asub and ',' not in Asub:
    #                    Asub = Asub + "-" + Asub

    #                if "(" in item[2] or "-" in item[2]:
    #                   Atext = '"' + item[2] + '"'
    #                else:
    #                    Atext = item[2]
    #                A_feature37 += '|' + item[0] + ':' + Asub + '(' + Atext + ')'
    #    else:
    #        A_feature37 = '-'
                
        # Get isoform info for A
        Aiso = get_isoforms(xint,conn)
        Atag = get_tag_info(xint,conn)

        if Aiso:
            A_annot26 = 'comment:"' + Aiso[0][0] + ' specific"'
            if Atag:
                 for each in Atag:
                    A_annot26 += '|comment:"' + each[0] + '"'
        elif len(Aiso) == 0:
            if Atag:
                A_annot26 = 'comment:"' + Atag[0][0] + '"'
                if len(Atag) > 1:
                    for each in Atag[1:]:
                        A_annot26 += '|comment:"' + each[0] + '"'
            else:
                A_annot26 = '-'

        # Get information for interactor B if there is one
        if len(parts) > 1:
            B_id2 = 'flybase:' + parts[1][0]
            if '(' in parts[1][1] or ':' in parts[1][1]:
                geneB = '"' + parts[1][1] + '"'
            else:
                geneB = parts[1][1]
            B_name6 = 'flybase:' + geneB + '(gene name)'
            if '-XP' in parts[1][3]:
                B_type22 = 'psi-mi:"MI:0326"(protein)'
            elif '-XR' in parts[1][3]:
                B_type22 = 'psi-mi:"MI:0320"(ribonucleic acid)'
            else:
                B_type22 = 'psi-mi:"MI:0329"(unknown participant)'

            gid = (parts[1][0],)

            CG = get_CG_id(gid,conn)
            if parts[1][0] in gene_map:
                Entrez = gene_map[parts[1][0]]
            else:
                Entrez = get_Entrez_id(gid,conn)

            if CG and Entrez:
                B_altids4 = 'flybase:' + CG[0][0] + '|entrez gene/locuslink:' + Entrez
            elif Entrez and not CG:
                B_altids4 = 'entrez gene/locuslink:' + Entrez
            elif CG and not Entrez:
                B_altids4 = 'flybase:' + CG[0][0]
            else:
                B_altids4 = '-'

            xid = (parts[1][2],)
            tax = get_taxid(xid,conn)
            if tax is not None:
                B_tax11 = 'taxid:' + tax[0] + '("' + tax[1] + ' ' + tax[2] + '")'
            else:
                 B_tax11 = '-'
            xint = (parts[1][2],int[0],)
            rol = get_role(xint,conn)
            if rol[0][2] in role_list:
                B_role20 =  'psi-mi:"MI:' + rol[0][0] + '"(' + rol[0][1] + ')'
                B_role18 = '-'
            else:
                B_role18 = 'psi-mi:"MI:' + rol[0][0] + '"(' + rol[0][1] + ')'
                B_role20 = 'psi-mi:"MI:0499"(unspecified role)'


    #    Comment out section to get subregions until data can be sufficiently standardized for mitab parser
    #    Get subregion/experimental feature info for B
    #        Bsubs = get_subregions(xint,conn)
    #        if Bsubs:
    #            if 'aa ' in Bsubs[0][1] or 'nt ' in Bsubs[0][1]:
    #                Bsubregion = Bsubs[0][1][3:]
    #                if 'aa ' in Bsubregion:
    #                    Bsubregion = Bsubregion.replace(" and aa ", ",")
    #            else:
    #                Bsubregion = Bsubs[0][1]
    #            if ('-' in Bsubregion and regexp.search(Bsubregion)) or ('-' not in Bsubregion and ',' in Bsubregion):
    #                components = Bsubregion.split(",")
    #                Bsubregion = components[0] + "-" + components[len(components)-1]
    #            elif '-' not in Bsubregion and ',' not in Bsubregion:
    #                Bsubregion = Bsubregion + "-" + Bsubregion
    #            if "(" in Bsubs[0][2] or "-" in Bsubs[0][2]:
    #                Btext = '"' + Bsubs[0][2] + '"'
    #            else:
    #                Btext = Bsubs[0][2]
    #            B_feature38 = Bsubs[0][0] + ':' + Bsubregion + '(' + Btext + ')'
    #            if len(Bsubs) > 1:
    #                for item in Bsubs[1:]:
    #                    if 'aa ' in item[1] or 'nt ' in item[1]:
    #                        Bsub = item[1][3:]
    #                        if 'aa ' in Bsub:
    #                            print(Bsub)
    #                            Bsub = Bsub.replace(" and aa ", ",")
    #                            print(Bsub)
    #                    else:
    #                        Bsub = item[1]
    #                    if ('-' in Bsub and regexp.search(Bsub)) or ('-' not in Bsub and ',' in Bsub):
    #                        components = Bsub.split(",")
    #                        Bsub = components[0] + "-" + components[len(components)-1]
    #                    elif '-' not in Bsub and ',' not in Bsub:
    #                        Bsub = Bsub + "-" + Bsub
    #                    if "(" in item[2] or "-" in item[2]:
    #                        Btext = '"' + item[2] + '"'
    #                    else:
    #                        Btext = item[2]
    #                    B_feature38 += '|' + item[0] + ':' + Bsub + '(' + Btext + ')'
    #        else:
    #            B_feature38 = '-'

            Biso = get_isoforms(xint,conn)
            Btag = get_tag_info(xint,conn)
            if Biso:
                B_annot27 = 'comment:"' + Biso[0][0] + ' specific"'
                if Btag:
                    for each in Btag:
                        B_annot27 += '|comment:"' + each[0] + '"'
            elif len(Biso) == 0:
                if len(Btag) >= 1:
                    B_annot27 = 'comment:"' + Btag[0][0] + '"'
                    if len(Btag) > 1:
                        for each in Btag[1:]:
                            B_annot27 += '|comment:"' + each[0] + '"'
                else:
                    B_annot27 = '-'

        # Otherwise fill B fields with A values
        else:
            B_id2 = A_id1
            B_name6 = A_name5
            B_altids4 = A_altids3
            B_tax11 = A_tax10
            B_type22 = A_type21
            B_role20 = A_role19
            B_role18 = A_role17
    #        B_feature38 = A_feature37
            B_annot27 = A_annot26
                
        # Print line for each interaction
        with open(filename, 'a') as csvfile:
            csv_writer = csv.writer(csvfile, quotechar = '', quoting=csv.QUOTE_NONE, delimiter = '\t')
            csv_writer.writerow([A_id1, B_id2, A_altids3, B_altids4,
                A_name5, B_name6, assays7, ref8, pubid9, A_tax10, B_tax11,
                intype12, 'psi-mi:"MI:0478"(flybase)', int_id14, '-',
                '-', A_role17, B_role18, A_role19, B_role20, A_type21,
                B_type22, '-', '-', int_xref25, A_annot26, B_annot27, annots28, '-',
                '-', '-', '-', '-', '-', '-', 'false', '-', '-', '-',
                '-', '-', '-'])


def main():
    parser = argparse.ArgumentParser(description='Query Chado.')
    parser.add_argument('-s', '--pgserver', help='Postgres server', required=True)
    parser.add_argument('-d', '--database', help='Reporting database', required=True)
    parser.add_argument('-u', '--username', help='Postgres User name', required=True)
    parser.add_argument('-p', '--password', help='Postgres password', required=True)
    parser.add_argument('-r', '--release', help='FlyBase release used', required=True)

    args = parser.parse_args() 
    server = args.pgserver
    database = args.database
    username = args.username
    password = args.password
    release = args.release
    
    # Define connection
    conn_string = "host='%s' dbname='%s' user='%s' password='%s'" % (server, database, username, password)
    filename = "physical_interactions_mitab_%s_reporting.tsv" % (release)

    # Attempt to get a connection
    conn = psycopg2.connect(conn_string)
    query_for_ints(filename,conn)

    # Close the connection
    conn.close()

if __name__ == "__main__":
    main()